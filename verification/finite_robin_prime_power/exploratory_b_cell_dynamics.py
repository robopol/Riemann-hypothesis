#!/usr/bin/env python3
"""Explore the prime-power cell dynamics of the Weil function B(t).

This is an exploratory binary64 computation.  It is not an interval
certificate, a proof of a polynomial envelope, or a proof of RH.

The script uses the exact finite prime-side formula for B(t), differentiates
it inside each open prime-power cell, verifies the derivative jump at every
event, locates cellwise stationary points, and tests several deliberately
simple Lyapunov/energy candidates.  The tests are intended to discover a
plausible invariant or to expose why a candidate cannot work.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import time
from pathlib import Path
from typing import Any

import numpy as np

from exploratory_integral_prime_side import (
    EULER_GAMMA,
    LOG_TWO_PI,
    WEIL_C,
    PrimeTables,
    build_prime_tables,
    prime_side_weil_b,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_REPORT = ROOT / "exploratory_b_cell_dynamics_report.json"


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def trivial_f_array(y: np.ndarray) -> np.ndarray:
    """Evaluate F(y)=2 atanh(y)-2y+y log(1-y^2) stably."""

    y = np.asarray(y, dtype=np.float64)
    output = np.empty_like(y)
    small = y < 0.01
    if np.any(small):
        ys = y[small]
        square = ys * ys
        power = ys * square
        total = np.zeros_like(ys)
        # Six terms make the omitted tail negligible for y < 0.01.
        for layer in range(1, 7):
            total -= power / (layer * (2 * layer + 1))
            power *= square
        output[small] = total
    if np.any(~small):
        yl = y[~small]
        output[~small] = (
            (1.0 + yl) * np.log1p(yl)
            + (yl - 1.0) * np.log1p(-yl)
            - 2.0 * yl
        )
    return output


def cell_b_and_derivative(
    t: np.ndarray | float,
    s_lambda: np.ndarray | float,
    psi: np.ndarray | float,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate B and B' with fixed cell constants S_Lambda and psi."""

    t_array = np.asarray(t, dtype=np.float64)
    s_array = np.asarray(s_lambda, dtype=np.float64)
    psi_array = np.asarray(psi, dtype=np.float64)
    root_x = np.exp(0.5 * t_array)
    y = np.exp(-t_array)
    f_value = trivial_f_array(y)
    f_prime = np.log1p(-(y * y))
    a_value = s_array + EULER_GAMMA + 1.0
    q_value = psi_array + LOG_TWO_PI

    b_value = root_x * (a_value - t_array - 0.5 * f_value) - q_value / root_x
    derivative = root_x * (
        0.5 * (a_value - t_array)
        - 1.0
        - 0.25 * f_value
        + 0.5 * y * f_prime
    ) + 0.5 * q_value / root_x
    return b_value, derivative


def scalar_cell_derivative(t: float, s_lambda: float, psi: float) -> float:
    """Scalar wrapper used by stationary-point bisection."""

    return float(cell_b_and_derivative(t, s_lambda, psi)[1])


def stationary_point(
    left: float,
    right: float,
    s_lambda: float,
    psi: float,
    derivative_left: float,
    derivative_right: float,
) -> tuple[float, float]:
    """Locate a derivative zero in a cell by safeguarded bisection."""

    if derivative_left == 0.0:
        root = left
    elif derivative_right == 0.0:
        root = right
    else:
        lo = left
        hi = right
        f_lo = derivative_left
        for _ in range(55):
            midpoint = 0.5 * (lo + hi)
            f_mid = scalar_cell_derivative(midpoint, s_lambda, psi)
            if f_mid == 0.0:
                lo = hi = midpoint
                break
            if (f_mid > 0.0) == (f_lo > 0.0):
                lo = midpoint
                f_lo = f_mid
            else:
                hi = midpoint
        root = 0.5 * (lo + hi)
    b_value = float(cell_b_and_derivative(root, s_lambda, psi)[0])
    return root, b_value


def extremum_row(
    t: float, b_value: float, derivative: float | None, kind: str, event: int | None
) -> dict[str, Any]:
    """Format a state/extremum row for JSON."""

    row: dict[str, Any] = {
        "kind": kind,
        "t": float(t),
        "x": float(math.exp(t)),
        "B": float(b_value),
        "abs_B": float(abs(b_value)),
    }
    if derivative is not None:
        row["B_prime"] = float(derivative)
    if event is not None:
        row["prime_power_event"] = int(event)
    return row


def summarize_values(
    t: np.ndarray,
    values: np.ndarray,
    events: np.ndarray,
    label: str,
) -> dict[str, Any]:
    """Return signed and absolute extrema for an event-valued array."""

    minimum = int(np.argmin(values))
    maximum = int(np.argmax(values))
    absolute = int(np.argmax(np.abs(values)))
    return {
        "label": label,
        "minimum": extremum_row(
            float(t[minimum]), float(values[minimum]), None, "event", int(events[minimum])
        ),
        "maximum": extremum_row(
            float(t[maximum]), float(values[maximum]), None, "event", int(events[maximum])
        ),
        "maximum_absolute": extremum_row(
            float(t[absolute]), float(values[absolute]), None, "event", int(events[absolute])
        ),
    }


def logarithmic_bin_diagnostics(
    t: np.ndarray,
    b_value: np.ndarray,
    derivative_left: np.ndarray,
    derivative_right: np.ndarray,
    bin_width: float,
) -> list[dict[str, Any]]:
    """Measure maxima in fixed-width logarithmic bins."""

    bin_index = np.floor(t / bin_width).astype(np.int64)
    rows: list[dict[str, Any]] = []
    for index in range(int(bin_index[0]), int(bin_index[-1]) + 1):
        mask = bin_index == index
        if not np.any(mask):
            continue
        local_t = t[mask]
        local_b = b_value[mask]
        local_d = np.maximum(np.abs(derivative_left[mask]), np.abs(derivative_right[mask]))
        rows.append(
            {
                "t_left": float(index * bin_width),
                "t_right": float((index + 1) * bin_width),
                "event_count": int(np.count_nonzero(mask)),
                "max_abs_B_at_events": float(np.max(np.abs(local_b))),
                "max_abs_B_prime_one_sided": float(np.max(local_d)),
                "rms_B_at_events": float(np.sqrt(np.mean(local_b * local_b))),
                "rms_B_prime_one_sided": float(np.sqrt(np.mean(local_d * local_d))),
                "last_t": float(local_t[-1]),
            }
        )
    return rows


def polynomial_envelope_tests(
    t: np.ndarray,
    b_value: np.ndarray,
    derivative_left: np.ndarray,
    derivative_right: np.ndarray,
    powers: list[float],
) -> list[dict[str, Any]]:
    """Test finite-scan polynomial normalizations for B and one-sided B'."""

    rows: list[dict[str, Any]] = []
    derivative_abs = np.maximum(np.abs(derivative_left), np.abs(derivative_right))
    for power in powers:
        scale = np.power(1.0 + t, power)
        normalized_b = np.abs(b_value) / scale
        normalized_d = derivative_abs / scale
        b_index = int(np.argmax(normalized_b))
        d_index = int(np.argmax(normalized_d))
        rows.append(
            {
                "power": power,
                "max_abs_B_over_one_plus_t_power": float(normalized_b[b_index]),
                "B_maximizer_t": float(t[b_index]),
                "B_maximizer_x": float(math.exp(float(t[b_index]))),
                "max_abs_B_prime_over_one_plus_t_power": float(normalized_d[d_index]),
                "B_prime_maximizer_t": float(t[d_index]),
                "B_prime_maximizer_x": float(math.exp(float(t[d_index]))),
            }
        )
    return rows


def energy_candidate_tests(
    t: np.ndarray,
    b_value: np.ndarray,
    derivative_left: np.ndarray,
    derivative_right: np.ndarray,
) -> list[dict[str, Any]]:
    """Test positive quadratic energies at jumps and across open cells."""

    specifications = [
        ("B2_plus_1e-4_D2", 1.0e-4, 0.0, 0.0),
        ("B2_plus_1e-3_D2", 1.0e-3, 0.0, 0.0),
        ("B2_plus_1e-2_D2", 1.0e-2, 0.0, 0.0),
        ("B2_plus_1e-1_D2", 1.0e-1, 0.0, 0.0),
        ("B2_plus_D2", 1.0, 0.0, 0.0),
        ("B2_plus_1e-2_D_plus_B2", 1.0e-2, 1.0, 0.0),
        ("B2_plus_1e-2_D_minus_B2", 1.0e-2, -1.0, 0.0),
        ("B2_plus_D2_over_one_plus_t2", 1.0, 0.0, 1.0),
        ("B2_plus_D2_over_one_plus_t4", 1.0, 0.0, 2.0),
    ]
    rows: list[dict[str, Any]] = []
    for name, coefficient, shift, scale_power in specifications:
        scale = np.power(1.0 + t, scale_power)
        shifted_left = (derivative_left + shift * b_value) / scale
        shifted_right = (derivative_right + shift * b_value) / scale
        energy_left = b_value * b_value + coefficient * shifted_left * shifted_left
        energy_right = b_value * b_value + coefficient * shifted_right * shifted_right
        jump_change = energy_right - energy_left

        # The right state at event i evolves continuously to the left state at i+1.
        cell_change = energy_left[1:] - energy_right[:-1]
        jump_positive = jump_change > 1e-15
        cell_positive = cell_change > 1e-15
        total_change = energy_right[1:] - energy_right[:-1]
        total_positive = total_change > 1e-15
        worst_jump = int(np.argmax(jump_change))
        worst_cell = int(np.argmax(cell_change))
        worst_total = int(np.argmax(total_change))
        max_energy = int(np.argmax(np.maximum(energy_left, energy_right)))
        rows.append(
            {
                "name": name,
                "definition": (
                    f"E=B^2+{coefficient:g}*((B'+{shift:g}*B)/(1+t)^{scale_power:g})^2"
                ),
                "event_jump_increase_count": int(np.count_nonzero(jump_positive)),
                "event_jump_increase_fraction": float(np.mean(jump_positive)),
                "open_cell_increase_count": int(np.count_nonzero(cell_positive)),
                "open_cell_increase_fraction": float(np.mean(cell_positive)),
                "event_to_event_increase_count": int(np.count_nonzero(total_positive)),
                "event_to_event_increase_fraction": float(np.mean(total_positive)),
                "largest_event_jump_increase": float(jump_change[worst_jump]),
                "largest_event_jump_increase_t": float(t[worst_jump]),
                "largest_open_cell_increase": float(cell_change[worst_cell]),
                "largest_open_cell_increase_from_t": float(t[worst_cell]),
                "largest_event_to_event_increase": float(total_change[worst_total]),
                "largest_event_to_event_increase_from_t": float(t[worst_total]),
                "maximum_energy": float(
                    max(energy_left[max_energy], energy_right[max_energy])
                ),
                "maximum_energy_t": float(t[max_energy]),
                "endpoint_energy_right": float(energy_right[-1]),
                "classification": (
                    "monotone_on_scanned_events_and_cells"
                    if not np.any(jump_positive) and not np.any(cell_positive)
                    else "fails_monotonicity_on_finite_scan"
                ),
            }
        )
    return rows


def constant_quadratic_jump_feasibility(
    t: np.ndarray,
    b_value: np.ndarray,
    derivative_left: np.ndarray,
    event_jump: np.ndarray,
    events: np.ndarray,
) -> dict[str, Any]:
    """Test every constant positive-definite quadratic form at event jumps.

    For E=u*B^2+2*v*B*B'+w*B'^2 with w>0, an event jump j>0 has
    Delta E/(w*j)=2*(v/w)*B+2*B'_-+j.  Thus all jump inequalities reduce to
    a one-dimensional feasibility interval for r=v/w.  The coefficient u can
    always be chosen larger than r^2*w to make the form positive definite.
    """

    threshold = 1e-14
    positive = b_value > threshold
    negative = b_value < -threshold
    near_zero = ~(positive | negative)
    boundary = np.full_like(b_value, np.nan)
    np.divide(
        -(2.0 * derivative_left + event_jump),
        2.0 * b_value,
        out=boundary,
        where=~near_zero,
    )

    upper_values = np.where(positive, boundary, np.inf)
    lower_values = np.where(negative, boundary, -np.inf)
    upper_index = int(np.argmin(upper_values))
    lower_index = int(np.argmax(lower_values))
    upper = float(upper_values[upper_index])
    lower = float(lower_values[lower_index])
    near_zero_violation = near_zero & ((2.0 * derivative_left + event_jump) > 0.0)
    feasible = lower <= upper and not np.any(near_zero_violation)

    return {
        "family": "E=u*B^2+2*v*B*B'+w*B'^2, w>0, u*w-v^2>0",
        "jump_condition": (
            "With r=v/w, every event requires 2*r*B+2*B'_-+j<=0."
        ),
        "required_r_lower_bound_from_negative_B_events": lower,
        "lower_bound_witness": {
            "event": int(events[lower_index]),
            "t": float(t[lower_index]),
            "B": float(b_value[lower_index]),
            "B_prime_left": float(derivative_left[lower_index]),
            "jump": float(event_jump[lower_index]),
        },
        "required_r_upper_bound_from_positive_B_events": upper,
        "upper_bound_witness": {
            "event": int(events[upper_index]),
            "t": float(t[upper_index]),
            "B": float(b_value[upper_index]),
            "B_prime_left": float(derivative_left[upper_index]),
            "jump": float(event_jump[upper_index]),
        },
        "feasibility_gap_upper_minus_lower": upper - lower,
        "near_zero_B_event_count": int(np.count_nonzero(near_zero)),
        "near_zero_constraint_violation_count": int(
            np.count_nonzero(near_zero_violation)
        ),
        "constant_positive_definite_quadratic_can_be_nonincreasing_at_all_scanned_jumps": bool(
            feasible
        ),
        "interpretation": (
            "An empty r interval rigorously disqualifies this entire constant "
            "quadratic family on the finite binary64 data (subject to numerical "
            "rounding); it does not disqualify time-dependent or prime-adapted energies."
        ),
    }


def regularized_bhat_energy_report(
    t: np.ndarray,
    s_left: np.ndarray,
    s_right: np.ndarray,
    psi_left: np.ndarray,
    psi_right: np.ndarray,
    b_left: np.ndarray,
    b_right: np.ndarray,
    event_jump: np.ndarray,
    events: np.ndarray,
) -> dict[str, Any]:
    """Test the exact near-degenerate polynomial energy family for Bhat."""

    x = np.exp(t)
    root_x = np.sqrt(x)
    trivial_correction = 0.5 * root_x * trivial_f_array(1.0 / x)

    def states(
        s_lambda: np.ndarray, psi: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        u_value = s_lambda - t + EULER_GAMMA
        v_value = psi - x + LOG_TWO_PI
        bhat = root_x * u_value - v_value / root_x
        bhat_prime = 0.5 * (root_x * u_value + v_value / root_x)
        return u_value, v_value, bhat, bhat_prime

    u_left, v_left, bhat_left, derivative_left = states(s_left, psi_left)
    u_right, v_right, bhat_right, derivative_right = states(s_right, psi_right)
    bhat_event = 0.5 * (bhat_left + bhat_right)
    correction_error_left = bhat_left - (b_left + trivial_correction)
    correction_error_right = bhat_right - (b_right + trivial_correction)
    derivative_jump_error = derivative_right - derivative_left - event_jump

    exponents = [0, 1, 2, 3, 4, 6, 8]
    rows: list[dict[str, Any]] = []
    tolerance = 2e-13
    for exponent in exponents:
        epsilon = np.power(1.0 + t, -float(exponent))

        def energy(derivative: np.ndarray) -> np.ndarray:
            return (
                (1.0 - 0.5 * epsilon) * bhat_event * bhat_event
                + 2.0 * epsilon * derivative * derivative
            )

        energy_left = energy(derivative_left)
        energy_right = energy(derivative_right)
        jump_change = energy_right - energy_left
        predicted_jump = (
            4.0 * epsilon * event_jump * derivative_left
            + 2.0 * epsilon * event_jump * event_jump
        )
        cell_change = energy_left[1:] - energy_right[:-1]
        event_to_event = energy_right[1:] - energy_right[:-1]
        jump_increase = predicted_jump > 0.0
        cell_increase = cell_change > 0.0
        event_to_event_increase = event_to_event > 0.0
        maximum_index = int(
            np.argmax(np.maximum(energy_left, energy_right))
        )
        normalized = np.maximum(energy_left, energy_right) / np.power(
            1.0 + t, float(exponent)
        )
        normalized_index = int(np.argmax(normalized))
        endpoint_normalized = energy_right[-1] / (1.0 + t[-1]) ** exponent
        initial_left_energy = float(energy_left[0])
        signed_jump_sum = float(np.sum(jump_change))
        signed_cell_sum = float(np.sum(cell_change))
        telescoping_error = (
            float(energy_right[-1] - energy_left[0])
            - signed_jump_sum
            - signed_cell_sum
        )
        rows.append(
            {
                "r": exponent,
                "epsilon": "epsilon(t)=(1+t)^(-r)",
                "energy": (
                    "E=(1-epsilon/2)*Bhat^2+2*epsilon*Bhat'^2"
                ),
                "equivalent_UV_form": (
                    "E=x*U^2+V^2/x-2*(1-epsilon)*U*V"
                ),
                "jump_increase_count": int(
                    np.count_nonzero(jump_increase)
                ),
                "jump_increase_fraction": float(
                    np.mean(jump_increase)
                ),
                "open_cell_increase_count": int(
                    np.count_nonzero(cell_increase)
                ),
                "open_cell_increase_fraction": float(
                    np.mean(cell_increase)
                ),
                "event_to_event_increase_count": int(
                    np.count_nonzero(event_to_event_increase)
                ),
                "event_to_event_increase_fraction": float(
                    np.mean(event_to_event_increase)
                ),
                "counting_note": (
                    "Primary increase counts use the raw binary64 sign; robust "
                    "counts below require an absolute increase greater than 2e-13."
                ),
                "robust_jump_increase_count_above_2e_minus_13": int(
                    np.count_nonzero(predicted_jump > tolerance)
                ),
                "robust_open_cell_increase_count_above_2e_minus_13": int(
                    np.count_nonzero(cell_change > tolerance)
                ),
                "robust_event_to_event_increase_count_above_2e_minus_13": int(
                    np.count_nonzero(event_to_event > tolerance)
                ),
                "sum_positive_jump_increments": float(
                    np.sum(np.maximum(jump_change, 0.0))
                ),
                "sum_negative_jump_increments": float(
                    np.sum(np.minimum(jump_change, 0.0))
                ),
                "signed_jump_sum": signed_jump_sum,
                "signed_open_cell_drift_sum": signed_cell_sum,
                "signed_jump_plus_cell_sum": signed_jump_sum + signed_cell_sum,
                "endpoint_minus_initial_left_energy": float(
                    energy_right[-1] - energy_left[0]
                ),
                "telescoping_balance_error": telescoping_error,
                "maximum_one_sided_event_energy": float(
                    max(energy_left[maximum_index], energy_right[maximum_index])
                ),
                "maximum_energy_t": float(t[maximum_index]),
                "maximum_energy_event": int(events[maximum_index]),
                "endpoint_right_energy": float(energy_right[-1]),
                "initial_left_energy_before_event_2": initial_left_energy,
                "maximum_energy_over_initial_left_energy": float(
                    max(energy_left[maximum_index], energy_right[maximum_index])
                    / initial_left_energy
                ),
                "endpoint_energy_over_initial_left_energy": float(
                    energy_right[-1] / initial_left_energy
                ),
                "maximum_energy_divided_by_one_plus_t_power_r": float(
                    normalized[normalized_index]
                ),
                "normalized_maximum_t": float(t[normalized_index]),
                "normalized_maximum_event": int(events[normalized_index]),
                "endpoint_energy_divided_by_one_plus_t_power_r": float(
                    endpoint_normalized
                ),
                "maximum_sqrt_energy_over_one_plus_t_power_r_over_2": float(
                    math.sqrt(float(normalized[normalized_index]))
                ),
                "max_abs_exact_jump_formula_error": float(
                    np.max(np.abs(jump_change - predicted_jump))
                ),
            }
        )

    return {
        "status": (
            "Finite binary64 tests of an exact energy family. Event counts and "
            "sums do not establish a uniform polynomial envelope."
        ),
        "regularization": {
            "T": "T(y)=F(y)=2*atanh(y)-2*y+y*log(1-y^2)",
            "Bhat_from_B": "Bhat(t)=B(t)+exp(t/2)*T(exp(-t))/2",
            "U": "U(t)=S_Lambda(exp(t))-t+gamma_E",
            "V": "V(t)=psi(exp(t))-exp(t)+log(2pi)",
            "Bhat_prime_side": "Bhat=exp(t/2)*U-exp(-t/2)*V",
            "Bhat_prime": "Bhat'=(exp(t/2)*U+exp(-t/2)*V)/2",
            "event_jump": "Delta Bhat'=Lambda(n)/sqrt(n)",
            "cell_energy_derivative": (
                "For epsilon=(1+t)^(-r), E'=2*Bhat*Bhat'"
                "-4*epsilon*exp(t/2)*Bhat'"
                "+r*epsilon*Bhat^2/(2*(1+t))"
                "-2*r*epsilon*Bhat'^2/(1+t)."
            ),
        },
        "identity_validation": {
            "max_abs_Bhat_direct_minus_B_plus_trivial_left": float(
                np.max(np.abs(correction_error_left))
            ),
            "max_abs_Bhat_direct_minus_B_plus_trivial_right": float(
                np.max(np.abs(correction_error_right))
            ),
            "max_abs_Bhat_right_minus_left": float(
                np.max(np.abs(bhat_right - bhat_left))
            ),
            "max_abs_derivative_jump_error": float(
                np.max(np.abs(derivative_jump_error))
            ),
        },
        "Bhat_event_extrema": {
            "minimum": float(np.min(bhat_event)),
            "maximum": float(np.max(bhat_event)),
            "maximum_absolute": float(np.max(np.abs(bhat_event))),
            "maximum_abs_Bhat_prime_one_sided": float(
                np.max(
                    np.maximum(np.abs(derivative_left), np.abs(derivative_right))
                )
            ),
        },
        "rows": rows,
    }


def forcing_g_and_derivative(t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the exact cell forcing g and its positive derivative."""

    y = np.exp(-t)
    root_x = np.exp(0.5 * t)
    denominator = 1.0 - y * y
    quotient = y**3 / denominator
    g_value = root_x * (1.0 - quotient)
    g_prime = root_x * (
        0.5 * (1.0 - quotient)
        + (3.0 * y**3 - y**5) / (denominator * denominator)
    )
    return g_value, g_prime


def first_crossing_energy_report(
    t: np.ndarray,
    b_value: np.ndarray,
    derivative_left: np.ndarray,
    derivative_right: np.ndarray,
    event_jump: np.ndarray,
    events: np.ndarray,
) -> dict[str, Any]:
    """Analyze Phi=C-B and its exact first-crossing energy H."""

    phi = WEIL_C - b_value
    phi_prime_left = -derivative_left
    phi_prime_right = -derivative_right
    g_value, g_prime = forcing_g_and_derivative(t)
    k_value = g_value - 0.25 * WEIL_C

    def hamiltonian(phi_prime: np.ndarray) -> np.ndarray:
        return (
            phi_prime * phi_prime
            - 0.25 * phi * phi
            - 2.0 * k_value * phi
        )

    h_left = hamiltonian(phi_prime_left)
    h_right = hamiltonian(phi_prime_right)
    cone_radius_squared = 0.25 * phi * phi + 2.0 * k_value * phi
    cone_radius = np.sqrt(np.maximum(cone_radius_squared, 0.0))
    cone_margin_left = cone_radius - np.abs(phi_prime_left)
    cone_margin_right = cone_radius - np.abs(phi_prime_right)
    cone_ratio_left = np.abs(phi_prime_left) / cone_radius
    cone_ratio_right = np.abs(phi_prime_right) / cone_radius
    observed_event_change = h_right - h_left
    predicted_event_change = event_jump * event_jump - 2.0 * event_jump * phi_prime_left
    event_change_error = observed_event_change - predicted_event_change

    # In an open cell H'=-2*k'(t)*Phi.  When Phi>=0 and k'>0, H decreases.
    open_cell_change = h_left[1:] - h_right[:-1]
    positive_tolerance = 2e-10
    open_cell_increase = open_cell_change > positive_tolerance
    event_increase = observed_event_change > positive_tolerance
    event_decrease = observed_event_change < -positive_tolerance
    endpoint_change = h_right[1:] - h_right[:-1]

    minimum_phi_index = int(np.argmin(phi))
    maximum_phi_index = int(np.argmax(phi))
    minimum_h_index = int(np.argmin(np.minimum(h_left, h_right)))
    maximum_h_index = int(np.argmax(np.maximum(h_left, h_right)))
    largest_impulse_index = int(np.argmax(observed_event_change))
    largest_cell_drop_index = int(np.argmin(open_cell_change))
    smallest_cone_left = int(np.argmin(cone_margin_left))
    smallest_cone_right = int(np.argmin(cone_margin_right))
    lower_cone_margin_right = cone_radius + phi_prime_right
    pre_kick_ratio = (cone_radius + phi_prime_left) / event_jump
    smallest_lower_cone_right = int(np.argmin(lower_cone_margin_right))
    smallest_pre_kick_ratio = int(np.argmin(pre_kick_ratio))
    largest_cone_ratio = int(
        np.argmax(np.maximum(cone_ratio_left, cone_ratio_right))
    )

    bin_index = np.floor(t).astype(np.int64)
    balance_rows: list[dict[str, Any]] = []
    for index in range(int(bin_index[0]), int(bin_index[-1]) + 1):
        mask = bin_index == index
        if not np.any(mask):
            continue
        positions = np.flatnonzero(mask)
        start = int(positions[0])
        stop = int(positions[-1]) + 1
        local_event_changes = observed_event_change[start:stop]
        local_cell_changes = open_cell_change[start : min(stop, open_cell_change.size)]
        balance_rows.append(
            {
                "t_left": float(index),
                "t_right": float(index + 1),
                "event_count": int(stop - start),
                "sum_positive_event_increments": float(
                    np.sum(np.maximum(local_event_changes, 0.0))
                ),
                "sum_negative_event_increments": float(
                    np.sum(np.minimum(local_event_changes, 0.0))
                ),
                "sum_open_cell_changes": float(np.sum(local_cell_changes)),
                "net_right_state_change_within_index_range": float(
                    h_right[stop - 1] - h_right[start]
                ),
                "min_Phi_at_events": float(np.min(phi[start:stop])),
                "max_H_one_sided": float(
                    np.max(np.maximum(h_left[start:stop], h_right[start:stop]))
                ),
            }
        )

    return {
        "exact_dynamics": {
            "Phi_definition": "Phi(t)=C-B(t), C=B(0)=2+gamma_E-log(4pi)",
            "cell_ode": "Phi''(t)-Phi(t)/4=k(t)",
            "g": "g(t)=exp(t/2)*(1-exp(-3t)/(1-exp(-2t)))",
            "k": "k(t)=g(t)-C/4",
            "energy": "H=Phi'^2-Phi^2/4-2*k(t)*Phi",
            "open_cell_law": "H'=-2*k'(t)*Phi(t)",
            "event_law": (
                "At n=p^m, Phi'_+=Phi'_- - j_n with j_n=Lambda(n)/sqrt(n), "
                "Delta H=j_n^2-2*j_n*Phi'_-"
            ),
            "first_crossing_use": (
                "At a hypothetical first downward crossing Phi=0, H=Phi'^2>=0. "
                "A proof maintaining H<0 up to that point would exclude the crossing."
            ),
            "discrete_reduction": (
                "While Phi>=0, open-cell evolution cannot increase H. Therefore, "
                "after an initial check, it is enough to prove at every event that "
                "Phi'_+=Phi'_- - j_n remains above the lower cone boundary -R, "
                "where R=sqrt(Phi^2/4+2*k*Phi)."
            ),
        },
        "forcing_signs_on_scan": {
            "minimum_g": float(np.min(g_value)),
            "minimum_k": float(np.min(k_value)),
            "minimum_k_prime": float(np.min(g_prime)),
            "k_positive_at_all_scanned_events": bool(np.all(k_value > 0.0)),
            "k_prime_positive_at_all_scanned_events": bool(np.all(g_prime > 0.0)),
        },
        "Phi_at_events": {
            "minimum": {
                "t": float(t[minimum_phi_index]),
                "x": float(events[minimum_phi_index]),
                "Phi": float(phi[minimum_phi_index]),
                "B": float(b_value[minimum_phi_index]),
            },
            "maximum": {
                "t": float(t[maximum_phi_index]),
                "x": float(events[maximum_phi_index]),
                "Phi": float(phi[maximum_phi_index]),
                "B": float(b_value[maximum_phi_index]),
            },
            "nonnegative_at_all_scanned_events": bool(np.all(phi >= -1e-9)),
        },
        "H_one_sided_extrema": {
            "minimum": {
                "t": float(t[minimum_h_index]),
                "x": float(events[minimum_h_index]),
                "H_left": float(h_left[minimum_h_index]),
                "H_right": float(h_right[minimum_h_index]),
            },
            "maximum": {
                "t": float(t[maximum_h_index]),
                "x": float(events[maximum_h_index]),
                "H_left": float(h_left[maximum_h_index]),
                "H_right": float(h_right[maximum_h_index]),
            },
            "endpoint_H_right": float(h_right[-1]),
            "H_right_negative_at_all_scanned_events": bool(np.all(h_right < 0.0)),
            "H_both_sides_negative_at_all_scanned_events": bool(
                np.all(h_right < 0.0) and np.all(h_left < 0.0)
            ),
        },
        "negative_H_cone_margins": {
            "equivalence": (
                "H<0 iff |Phi'|<R, where R=sqrt(Phi^2/4+2*k*Phi), "
                "provided Phi>=0."
            ),
            "minimum_left_margin_R_minus_abs_Phi_prime": float(
                cone_margin_left[smallest_cone_left]
            ),
            "minimum_left_margin_t": float(t[smallest_cone_left]),
            "minimum_left_margin_event": int(events[smallest_cone_left]),
            "minimum_right_margin_R_minus_abs_Phi_prime": float(
                cone_margin_right[smallest_cone_right]
            ),
            "minimum_right_margin_t": float(t[smallest_cone_right]),
            "minimum_right_margin_event": int(events[smallest_cone_right]),
            "maximum_abs_Phi_prime_over_R": float(
                max(
                    cone_ratio_left[largest_cone_ratio],
                    cone_ratio_right[largest_cone_ratio],
                )
            ),
            "maximum_ratio_t": float(t[largest_cone_ratio]),
            "maximum_ratio_event": int(events[largest_cone_ratio]),
            "both_one_sided_states_strictly_inside_cone": bool(
                np.all(cone_margin_left > 0.0) and np.all(cone_margin_right > 0.0)
            ),
            "minimum_event_lower_boundary_margin_R_plus_Phi_prime_right": float(
                lower_cone_margin_right[smallest_lower_cone_right]
            ),
            "minimum_event_lower_boundary_margin_t": float(
                t[smallest_lower_cone_right]
            ),
            "minimum_event_lower_boundary_margin_event": int(
                events[smallest_lower_cone_right]
            ),
            "minimum_pre_kick_ratio_R_plus_Phi_prime_left_over_j": float(
                pre_kick_ratio[smallest_pre_kick_ratio]
            ),
            "minimum_pre_kick_ratio_t": float(t[smallest_pre_kick_ratio]),
            "minimum_pre_kick_ratio_event": int(events[smallest_pre_kick_ratio]),
            "universal_kick_target": (
                "Prove (R+Phi'_-)/j_n>1 at every prime-power event; this is "
                "equivalent to the post-kick lower-cone margin R+Phi'_+>0."
            ),
        },
        "event_impulses": {
            "increase_count": int(np.count_nonzero(event_increase)),
            "increase_fraction": float(np.mean(event_increase)),
            "decrease_count": int(np.count_nonzero(event_decrease)),
            "largest_increase": float(observed_event_change[largest_impulse_index]),
            "largest_increase_t": float(t[largest_impulse_index]),
            "largest_increase_event": int(events[largest_impulse_index]),
            "sum_positive_increments": float(
                np.sum(np.maximum(observed_event_change, 0.0))
            ),
            "sum_negative_increments": float(
                np.sum(np.minimum(observed_event_change, 0.0))
            ),
            "max_abs_exact_jump_identity_error": float(
                np.max(np.abs(event_change_error))
            ),
        },
        "open_cells": {
            "increase_count_above_2e_minus_10": int(
                np.count_nonzero(open_cell_increase)
            ),
            "increase_fraction": float(np.mean(open_cell_increase)),
            "largest_apparent_increase": float(np.max(open_cell_change)),
            "largest_drop": float(open_cell_change[largest_cell_drop_index]),
            "largest_drop_from_t": float(t[largest_cell_drop_index]),
            "sum_changes": float(np.sum(open_cell_change)),
            "sign_prediction_holds_at_event_endpoints": bool(
                not np.any(open_cell_increase)
            ),
        },
        "event_to_event_H": {
            "increase_count_above_2e_minus_10": int(
                np.count_nonzero(endpoint_change > positive_tolerance)
            ),
            "increase_fraction": float(
                np.mean(endpoint_change > positive_tolerance)
            ),
            "largest_increase": float(np.max(endpoint_change)),
            "largest_decrease": float(np.min(endpoint_change)),
        },
        "log_bin_impulse_drift_balance": balance_rows,
    }


def locate_stationary_extrema(
    t: np.ndarray,
    s_lambda: np.ndarray,
    psi: np.ndarray,
    b_event: np.ndarray,
    derivative_left: np.ndarray,
    derivative_right: np.ndarray,
    events: np.ndarray,
) -> dict[str, Any]:
    """Locate derivative sign changes inside all scanned cells."""

    signs_cross = (
        (derivative_right[:-1] == 0.0)
        | (derivative_left[1:] == 0.0)
        | (np.signbit(derivative_right[:-1]) != np.signbit(derivative_left[1:]))
    )
    candidates = np.flatnonzero(signs_cross)
    stationary_rows: list[dict[str, Any]] = []
    minimum = extremum_row(
        float(t[int(np.argmin(b_event))]),
        float(np.min(b_event)),
        None,
        "prime_power_event",
        int(events[int(np.argmin(b_event))]),
    )
    maximum = extremum_row(
        float(t[int(np.argmax(b_event))]),
        float(np.max(b_event)),
        None,
        "prime_power_event",
        int(events[int(np.argmax(b_event))]),
    )
    max_abs = maximum if abs(maximum["B"]) >= abs(minimum["B"]) else minimum

    roots = np.empty(0, dtype=np.float64)
    values = np.empty(0, dtype=np.float64)
    if candidates.size:
        lo = t[candidates].copy()
        hi = t[candidates + 1].copy()
        f_lo = derivative_right[candidates].copy()
        cell_s = s_lambda[candidates]
        cell_psi = psi[candidates]
        # Vectorized bisection makes a complete multi-million-cell scan practical.
        for _ in range(48):
            midpoint = 0.5 * (lo + hi)
            f_mid = cell_b_and_derivative(midpoint, cell_s, cell_psi)[1]
            same_sign = np.signbit(f_mid) == np.signbit(f_lo)
            lo = np.where(same_sign, midpoint, lo)
            f_lo = np.where(same_sign, f_mid, f_lo)
            hi = np.where(same_sign, hi, midpoint)
        roots = 0.5 * (lo + hi)
        values = cell_b_and_derivative(roots, cell_s, cell_psi)[0]

        for root, value in zip(roots[:12], values[:12], strict=True):
            stationary_rows.append(
                extremum_row(
                    float(root), float(value), 0.0, "stationary_inside_cell", None
                )
            )

        stationary_minimum = int(np.argmin(values))
        stationary_maximum = int(np.argmax(values))
        stationary_absolute = int(np.argmax(np.abs(values)))
        candidate_rows = [
            extremum_row(
                float(roots[stationary_minimum]),
                float(values[stationary_minimum]),
                0.0,
                "stationary_inside_cell",
                None,
            ),
            extremum_row(
                float(roots[stationary_maximum]),
                float(values[stationary_maximum]),
                0.0,
                "stationary_inside_cell",
                None,
            ),
            extremum_row(
                float(roots[stationary_absolute]),
                float(values[stationary_absolute]),
                0.0,
                "stationary_inside_cell",
                None,
            ),
        ]
        for row in candidate_rows:
            if row["B"] < minimum["B"]:
                minimum = row
            if row["B"] > maximum["B"]:
                maximum = row
            if row["abs_B"] > max_abs["abs_B"]:
                max_abs = row

    return {
        "stationary_cell_count": int(candidates.size),
        "stationary_derivative_direction_counts": {
            "positive_to_negative": int(
                np.count_nonzero(
                    (derivative_right[candidates] > 0.0)
                    & (derivative_left[candidates + 1] < 0.0)
                )
            ),
            "negative_to_positive": int(
                np.count_nonzero(
                    (derivative_right[candidates] < 0.0)
                    & (derivative_left[candidates + 1] > 0.0)
                )
            ),
        },
        "first_stationary_points": stationary_rows,
        "global_minimum_over_events_and_stationary_points": minimum,
        "global_maximum_over_events_and_stationary_points": maximum,
        "global_maximum_absolute_over_events_and_stationary_points": max_abs,
        "all_extrema_respect_abs_B_le_B0_in_scan": bool(max_abs["abs_B"] <= WEIL_C + 1e-9),
    }


def sampled_direct_validation(
    tables: PrimeTables,
    event_t: np.ndarray,
    event_b: np.ndarray,
) -> dict[str, Any]:
    """Compare the vector event calculation with prime_side_weil_b."""

    indices = np.unique(
        np.linspace(0, event_t.size - 1, min(257, event_t.size), dtype=np.int64)
    )
    differences = []
    rows = []
    for index in indices:
        direct = prime_side_weil_b(tables, float(event_t[index]), "right")
        difference = float(event_b[index] - direct)
        differences.append(difference)
        if len(rows) < 10 or index == indices[-1]:
            rows.append(
                {
                    "event_index": int(index),
                    "prime_power": int(tables.event_value[index]),
                    "t": float(event_t[index]),
                    "cell_formula_B": float(event_b[index]),
                    "existing_prime_side_B": float(direct),
                    "difference": difference,
                }
            )
    return {
        "sample_count": int(indices.size),
        "max_abs_difference": float(np.max(np.abs(differences))),
        "selected_rows": rows,
    }


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    """Build the complete exploratory dynamics report."""

    started = time.perf_counter()
    tables = build_prime_tables(args.limit)
    event_mask = (
        (tables.event_value >= args.scan_min) & (tables.event_value <= args.limit)
    )
    all_indices = np.flatnonzero(event_mask)
    if all_indices.size < 2:
        raise ValueError("the scan requires at least two prime-power events")
    first = int(all_indices[0])
    last = int(all_indices[-1]) + 1

    events = tables.event_value[first:last]
    event_log_prime = tables.event_log_prime[first:last]
    t = np.log(events.astype(np.float64))
    s_right = tables.s_lambda_prefix[first:last]
    psi_right = tables.psi_prefix[first:last]
    s_left = np.empty_like(s_right)
    psi_left = np.empty_like(psi_right)
    if first:
        s_left[0] = tables.s_lambda_prefix[first - 1]
        psi_left[0] = tables.psi_prefix[first - 1]
    else:
        s_left[0] = 0.0
        psi_left[0] = 0.0
    s_left[1:] = s_right[:-1]
    psi_left[1:] = psi_right[:-1]

    b_left, derivative_left = cell_b_and_derivative(t, s_left, psi_left)
    b_right, derivative_right = cell_b_and_derivative(t, s_right, psi_right)
    b_event = 0.5 * (b_left + b_right)
    expected_jump = event_log_prime / np.sqrt(events.astype(np.float64))
    observed_jump = derivative_right - derivative_left

    stationary = locate_stationary_extrema(
        t,
        s_right,
        psi_right,
        b_event,
        derivative_left,
        derivative_right,
        events,
    )
    energy = energy_candidate_tests(t, b_event, derivative_left, derivative_right)
    quadratic_feasibility = constant_quadratic_jump_feasibility(
        t, b_event, derivative_left, expected_jump, events
    )
    bhat_energy = regularized_bhat_energy_report(
        t,
        s_left,
        s_right,
        psi_left,
        psi_right,
        b_left,
        b_right,
        expected_jump,
        events,
    )
    first_crossing = first_crossing_energy_report(
        t,
        b_event,
        derivative_left,
        derivative_right,
        expected_jump,
        events,
    )
    bins = logarithmic_bin_diagnostics(
        t, b_event, derivative_left, derivative_right, args.bin_width
    )
    envelope = polynomial_envelope_tests(
        t, b_event, derivative_left, derivative_right, args.powers
    )
    direct = sampled_direct_validation(tables, t, b_right)

    continuity_error = b_right - b_left
    jump_error = observed_jump - expected_jump
    checks = {
        "event_continuity_max_error_below_1e_minus_8": bool(
            np.max(np.abs(continuity_error)) < 1e-8
        ),
        "derivative_jump_max_error_below_1e_minus_8": bool(
            np.max(np.abs(jump_error)) < 1e-8
        ),
        "direct_prime_side_sample_max_error_below_1e_minus_10": bool(
            direct["max_abs_difference"] < 1e-10
        ),
        "stationary_extrema_respect_known_finite_abs_bound": stationary[
            "all_extrema_respect_abs_B_le_B0_in_scan"
        ],
        "first_crossing_event_jump_identity_below_1e_minus_8": bool(
            first_crossing["event_impulses"][
                "max_abs_exact_jump_identity_error"
            ]
            < 1e-8
        ),
        "Bhat_regularization_identity_below_1e_minus_8": bool(
            max(
                bhat_energy["identity_validation"][
                    "max_abs_Bhat_direct_minus_B_plus_trivial_left"
                ],
                bhat_energy["identity_validation"][
                    "max_abs_Bhat_direct_minus_B_plus_trivial_right"
                ],
                bhat_energy["identity_validation"][
                    "max_abs_Bhat_right_minus_left"
                ],
                bhat_energy["identity_validation"][
                    "max_abs_derivative_jump_error"
                ],
            )
            < 1e-8
        ),
        "Bhat_polynomial_energy_jump_formulas_below_1e_minus_8": bool(
            max(
                row["max_abs_exact_jump_formula_error"]
                for row in bhat_energy["rows"]
            )
            < 1e-8
        ),
    }
    checks["all_pass"] = all(checks.values())

    return {
        "status": (
            "Exploratory binary64 event scan. Exact formulas are distinguished "
            "from finite numerical evidence. This is not an interval certificate, "
            "a polynomial-envelope proof, or a proof of RH."
        ),
        "input": {
            "scan_min": args.scan_min,
            "scan_limit": args.limit,
            "prime_count": int(tables.primes.size),
            "all_prime_power_event_count": int(tables.event_value.size),
            "scanned_event_count": int(events.size),
            "first_event": int(events[0]),
            "last_event": int(events[-1]),
            "t_min": float(t[0]),
            "t_max": float(t[-1]),
            "bin_width": args.bin_width,
        },
        "exact_cell_dynamics": {
            "cell_constants": (
                "On a cell with fixed S=sum_{n<=x}Lambda(n)/n and P=psi(x), "
                "A=S+gamma_E+1, Q=P+log(2pi)."
            ),
            "B_formula": (
                "B(t)=exp(t/2)*(A-t-F(exp(-t))/2)-Q*exp(-t/2)"
            ),
            "B_prime_formula": (
                "B'(t)=exp(t/2)*((A-t)/2-1-F(y)/4+y*log(1-y^2)/2)"
                "+Q*exp(-t/2)/2, y=exp(-t)"
            ),
            "cell_ode": (
                "B''(t)-B(t)/4=-g(t), where "
                "g(t)=exp(t/2)*(1-exp(-3t)/(1-exp(-2t)))"
            ),
            "event_law": (
                "At n=p^m, B(log n+) = B(log n-) and "
                "B'(log n+)-B'(log n-)=Lambda(n)/sqrt(n)."
            ),
            "derivation_note": "F'(y)=log(1-y^2).",
            "first_crossing_concavity": (
                "If B(t)<=C at a putative first-crossing stage, then "
                "B''<=C/4-g=-k<0 in every open cell. Thus interior stationary "
                "points are maxima; prime-power upward jumps of B' can create "
                "event minima but not event maxima."
            ),
        },
        "event_law_validation": {
            "max_abs_B_right_minus_left": float(np.max(np.abs(continuity_error))),
            "max_abs_observed_minus_expected_derivative_jump": float(
                np.max(np.abs(jump_error))
            ),
            "largest_expected_derivative_jump": float(np.max(expected_jump)),
            "smallest_expected_derivative_jump": float(np.min(expected_jump)),
        },
        "existing_prime_side_validation": direct,
        "event_extrema": {
            "B_at_events": summarize_values(t, b_event, events, "B"),
            "B_prime_left": summarize_values(
                t, derivative_left, events, "B_prime_left"
            ),
            "B_prime_right": summarize_values(
                t, derivative_right, events, "B_prime_right"
            ),
        },
        "stationary_extrema": stationary,
        "Phi_first_crossing_energy": first_crossing,
        "polynomial_normalizations": envelope,
        "logarithmic_bin_diagnostics": bins,
        "quadratic_energy_candidates": {
            "interpretation": (
                "A candidate useful for a direct Lyapunov proof would need a "
                "rigorous global inequality. Failure on one scanned transition "
                "disqualifies that simple monotonicity claim; finite success would "
                "still not prove it."
            ),
            "rows": energy,
            "monotone_candidate_count": sum(
                row["classification"] == "monotone_on_scanned_events_and_cells"
                for row in energy
            ),
        },
        "constant_quadratic_jump_feasibility": quadratic_feasibility,
        "regularized_Bhat_polynomial_energy": bhat_energy,
        "checks": checks,
        "runtime": {
            "seconds": time.perf_counter() - started,
            "python": platform.python_version(),
            "numpy": np.__version__,
            "platform": platform.platform(),
            "script_sha256": sha256(Path(__file__).resolve()),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=56_048_351)
    parser.add_argument("--scan-min", type=int, default=2)
    parser.add_argument("--bin-width", type=float, default=1.0)
    parser.add_argument(
        "--powers", type=float, nargs="+", default=[0.0, 0.5, 1.0, 2.0, 4.0]
    )
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()
    if args.limit < 100 or args.scan_min < 2 or args.scan_min >= args.limit:
        raise ValueError("require 2 <= scan-min < limit and limit >= 100")
    if args.bin_width <= 0.0:
        raise ValueError("bin-width must be positive")

    report = build_report(args)
    if not report["checks"]["all_pass"]:
        raise AssertionError(f"exploratory consistency checks failed: {report['checks']}")
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"wrote {args.report}")
    print(f"all checks pass: {report['checks']['all_pass']}")
    print(
        "scanned events: "
        f"{report['input']['scanned_event_count']:,} through {report['input']['last_event']:,}"
    )
    print(
        "stationary cells: "
        f"{report['stationary_extrema']['stationary_cell_count']:,}"
    )
    print(
        "max |B| over events/stationary points: "
        f"{report['stationary_extrema']['global_maximum_absolute_over_events_and_stationary_points']['abs_B']:.12g}"
    )
    print(
        "simple monotone energy candidates: "
        f"{report['quadratic_energy_candidates']['monotone_candidate_count']}"
    )


if __name__ == "__main__":
    main()
