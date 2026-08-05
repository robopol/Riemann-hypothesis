#!/usr/bin/env python3
"""Scan the event-centered signed triangular finite-difference branch.

At every prime-power event q=p^a and for h=lambda*log(q)/sqrt(q), evaluate

    D_h = (Phi(t+h)-Phi(t-h))/(2*h),
    E_h = p_bar-D_h,
    U   = 4*Phi(t)*Phi(2*h)-(Phi(t+h)-Phi(t-h))^2,

and the RH-conditional Turan reserve

    M_T = R-sqrt(Phi(t)*Phi(2*h))/h+E_h-J/2.

The computation uses exact finite prime-side prefix formulas evaluated in
binary64.  It is exploratory: it is neither interval arithmetic nor a proof
of the Turan inequality, RH, or the event reserve inequality.
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

from exploratory_eventwise_reserve_scan import (
    EULER_GAMMA,
    LOG_TWO_PI,
    WEIL_C,
    higher_prime_powers,
    prime_array_in_segment,
    simple_primes,
    trivial_f_array,
)


SCRIPT_PATH = Path(__file__)
HELPER_PATH = SCRIPT_PATH.with_name("exploratory_eventwise_reserve_scan.py")
DEFAULT_REPORT = SCRIPT_PATH.with_name(
    "exploratory_signed_triangular_scan_report.json"
)
DEFAULT_LAMBDAS = (0.25, 0.30, 0.50, 1.00)


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def lambda_key(value: float) -> str:
    """Return a stable key for one window multiplier."""

    return f"lambda_{value:.12g}".replace(".", "p")


def required_support_limit(center_limit: int, maximum_lambda: float) -> int:
    """Return a safe event cutoff for all endpoint Phi evaluations.

    For h(t)=lambda*t*exp(-t/2), h is maximized at t=2 whenever that
    point lies in the center range.  The right endpoint q*exp(h(q)) is
    monotone for lambda<=exp(2).  For larger multipliers, the global h
    maximum gives a conservative bound for that endpoint as well.
    """

    t_left = math.log(2.0)
    t_right = math.log(float(center_limit))
    peak_t = min(max(2.0, t_left), t_right)
    maximum_h = maximum_lambda * peak_t * math.exp(-0.5 * peak_t)
    endpoint_h = maximum_lambda * t_right * math.exp(-0.5 * t_right)
    if maximum_lambda <= math.exp(2.0):
        maximum_x_plus = center_limit * math.exp(endpoint_h)
    else:
        maximum_x_plus = center_limit * math.exp(maximum_h)
    maximum_x_2h = math.exp(2.0 * maximum_h)
    return math.ceil(max(maximum_x_plus, maximum_x_2h)) + 2


def empty_sign_counts() -> dict[str, int]:
    """Create an empty sign counter."""

    return {"negative": 0, "zero_binary64": 0, "positive": 0}


def add_sign_counts(target: dict[str, int], values: np.ndarray) -> None:
    """Accumulate binary64 signs."""

    target["negative"] += int(np.count_nonzero(values < 0.0))
    target["zero_binary64"] += int(np.count_nonzero(values == 0.0))
    target["positive"] += int(np.count_nonzero(values > 0.0))


def update_extreme(
    target: dict[str, Any],
    name: str,
    value: float,
    row: dict[str, Any],
    *,
    maximum: bool = False,
) -> None:
    """Update a minimum or maximum witness."""

    current = target.get(name)
    better = current is None or (
        value > float(current["value"])
        if maximum
        else value < float(current["value"])
    )
    if better:
        target[name] = {"value": float(value), "row": row}


def build_event_support(
    support_limit: int, segment_span: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build the ordered prime-power event support through support_limit."""

    base_primes = simple_primes(math.isqrt(support_limit))
    prime_pieces: list[np.ndarray] = []
    low = 2
    while low <= support_limit:
        high = min(support_limit, low + segment_span - 1)
        prime_pieces.append(prime_array_in_segment(low, high, base_primes))
        low = high + 1
    primes = np.concatenate(prime_pieces)
    higher_q, higher_prime, higher_exponent, higher_log = higher_prime_powers(
        support_limit, base_primes
    )
    prime_logs = np.log(primes.astype(np.float64))
    event_q = np.concatenate((primes, higher_q))
    event_prime = np.concatenate((primes, higher_prime))
    event_exponent = np.concatenate(
        (
            np.ones(primes.size, dtype=np.int16),
            higher_exponent,
        )
    )
    event_lambda = np.concatenate((prime_logs, higher_log))
    order = np.argsort(event_q, kind="stable")
    return (
        event_q[order],
        event_prime[order],
        event_exponent[order],
        event_lambda[order],
    )


def phi_at_x(
    x: np.ndarray,
    event_q: np.ndarray,
    cumulative_psi: np.ndarray,
    cumulative_s_lambda: np.ndarray,
) -> np.ndarray:
    """Evaluate Phi(log(x)) from inclusive prime-power prefixes."""

    positions = np.searchsorted(event_q, x, side="right") - 1
    psi = np.zeros_like(x)
    s_lambda = np.zeros_like(x)
    present = positions >= 0
    if np.any(present):
        selected = positions[present]
        psi[present] = cumulative_psi[selected]
        s_lambda[present] = cumulative_s_lambda[selected]

    log_x = np.log(x)
    root_x = np.sqrt(x)
    y = 1.0 / x
    f_value = trivial_f_array(y)
    f_prime = np.log1p(-(y * y))
    u = (
        s_lambda
        - log_x
        + EULER_GAMMA
        + 0.5 * y * f_prime
        - 0.5 * f_value
    )
    v = psi - x + LOG_TWO_PI + 0.5 * f_prime
    coordinate_b = root_x * u - v / root_x
    return WEIL_C - coordinate_b


def compact_row(
    index: int,
    q: np.ndarray,
    prime: np.ndarray,
    exponent: np.ndarray,
    lambda_value: float,
    h: np.ndarray,
    phi: np.ndarray,
    phi_minus: np.ndarray,
    phi_plus: np.ndarray,
    phi_2h: np.ndarray,
    p_bar: np.ndarray,
    d_h: np.ndarray,
    e_h: np.ndarray,
    radius: np.ndarray,
    jump: np.ndarray,
    turan_bound: np.ndarray,
    turan_defect: np.ndarray,
    reserve: np.ndarray,
    actual_post_reserve: np.ndarray,
) -> dict[str, Any]:
    """Return a compact witness row."""

    log_q = math.log(int(q[index]))
    eta = log_q / math.sqrt(int(q[index]))
    discrepancy = float(e_h[index])
    bound = float(turan_bound[index])
    return {
        "q": int(q[index]),
        "prime": int(prime[index]),
        "exponent": int(exponent[index]),
        "lambda": lambda_value,
        "t": log_q,
        "h": float(h[index]),
        "eta": eta,
        "Phi_t": float(phi[index]),
        "Phi_t_minus_h": float(phi_minus[index]),
        "Phi_t_plus_h": float(phi_plus[index]),
        "Phi_2h": float(phi_2h[index]),
        "p_bar": float(p_bar[index]),
        "D_h": float(d_h[index]),
        "E_h": discrepancy,
        "J": float(jump[index]),
        "R": float(radius[index]),
        "turan_secant_bound": bound,
        "turan_defect_U": float(turan_defect[index]),
        "conditional_reserve_M_T": float(reserve[index]),
        "actual_post_reserve": float(actual_post_reserve[index]),
        "conditional_gap_actual_minus_M_T": float(
            actual_post_reserve[index] - reserve[index]
        ),
        "E_h_over_J": discrepancy / float(jump[index]),
        "E_h_over_eta": discrepancy / eta,
        "E_h_over_R": discrepancy / float(radius[index]),
        "E_h_over_turan_bound": discrepancy / bound,
        "triangular_T_equal_2hE": float(2.0 * h[index] * e_h[index]),
    }


def new_decade(decade: int) -> dict[str, Any]:
    """Create a decade aggregation record."""

    return {
        "decade": decade,
        "q_left": 2 if decade == 0 else 10**decade,
        "q_right": 10 ** (decade + 1) - 1,
        "event_count": 0,
        "E_h_sign_counts": empty_sign_counts(),
        "M_T_sign_counts": empty_sign_counts(),
        "raw_negative_U_count": 0,
        "U_below_minus_1e_minus_12_count": 0,
        "extrema": {},
    }


def new_lambda_result(lambda_value: float) -> dict[str, Any]:
    """Create online aggregation state for one lambda."""

    return {
        "lambda": lambda_value,
        "event_count": 0,
        "E_h_sign_counts": empty_sign_counts(),
        "M_T_sign_counts": empty_sign_counts(),
        "raw_negative_U_count": 0,
        "U_below_minus_1e_minus_12_count": 0,
        "nonpositive_phi_2h_count": 0,
        "nonpositive_radius_squared_count": 0,
        "global_extrema": {},
        "decades": {},
        "maximum_algebra_identity_residual": 0.0,
        "maximum_post_slope_identity_residual": 0.0,
    }


def scan(args: argparse.Namespace) -> dict[str, Any]:
    """Run the signed triangular scan."""

    started = time.perf_counter()
    lambda_values = tuple(float(value) for value in args.lambdas)
    if any(value <= 0.0 or value >= math.sqrt(2.0) for value in lambda_values):
        raise ValueError(
            "lambda values must satisfy 0<lambda<sqrt(2) so every window has 0<h<t"
        )
    maximum_lambda = max(lambda_values)
    support_limit = required_support_limit(args.limit, maximum_lambda)
    event_q, event_prime, event_exponent, event_lambda = build_event_support(
        support_limit, args.segment_span
    )
    cumulative_psi = np.cumsum(event_lambda, dtype=np.float64)
    cumulative_s_lambda = np.cumsum(
        event_lambda / event_q.astype(np.float64), dtype=np.float64
    )
    center_count = int(np.searchsorted(event_q, args.limit, side="right"))

    lambda_results = {
        lambda_key(value): new_lambda_result(value) for value in lambda_values
    }
    processed_chunks = 0
    for start in range(0, center_count, args.chunk_size):
        stop = min(center_count, start + args.chunk_size)
        q = event_q[start:stop]
        prime = event_prime[start:stop]
        exponent = event_exponent[start:stop]
        lambda_q = event_lambda[start:stop]
        x = q.astype(np.float64)
        root_x = np.sqrt(x)
        log_x = np.log(x)
        y = 1.0 / x

        psi_right = cumulative_psi[start:stop]
        s_right = cumulative_s_lambda[start:stop]
        f_value = trivial_f_array(y)
        f_prime = np.log1p(-(y * y))
        u_right = (
            s_right
            - log_x
            + EULER_GAMMA
            + 0.5 * y * f_prime
            - 0.5 * f_value
        )
        v_right = psi_right - x + LOG_TWO_PI + 0.5 * f_prime
        coordinate_b = root_x * u_right - v_right / root_x
        phi = WEIL_C - coordinate_b
        p_right = -0.5 * (root_x * u_right + v_right / root_x)
        jump = lambda_q / root_x
        p_bar = p_right + 0.5 * jump
        b_factor = 2.0 * (1.0 - 1.0 / (x * (x * x - 1.0)))
        radius_squared = phi * (
            root_x * b_factor - 0.5 * WEIL_C + 0.25 * phi
        )
        radius = np.sqrt(np.maximum(radius_squared, 0.0))
        actual_post_reserve = radius + p_right
        decade_keys = np.floor(np.log10(x)).astype(np.int16)

        for lambda_value in lambda_values:
            key = lambda_key(lambda_value)
            result = lambda_results[key]
            h = lambda_value * log_x / root_x
            x_minus = x * np.exp(-h)
            x_plus = x * np.exp(h)
            x_2h = np.exp(2.0 * h)
            phi_minus = phi_at_x(
                x_minus, event_q, cumulative_psi, cumulative_s_lambda
            )
            phi_plus = phi_at_x(
                x_plus, event_q, cumulative_psi, cumulative_s_lambda
            )
            phi_2h = phi_at_x(
                x_2h, event_q, cumulative_psi, cumulative_s_lambda
            )
            difference = phi_plus - phi_minus
            d_h = difference / (2.0 * h)
            e_h = p_bar - d_h
            phi_product = phi * phi_2h
            turan_bound = np.sqrt(np.maximum(phi_product, 0.0)) / h
            turan_defect = 4.0 * phi_product - difference * difference
            reserve = radius - turan_bound + e_h - 0.5 * jump
            conditional_gap = actual_post_reserve - reserve
            algebra_residual = conditional_gap - (d_h + turan_bound)
            result["maximum_algebra_identity_residual"] = max(
                result["maximum_algebra_identity_residual"],
                float(np.max(np.abs(algebra_residual))),
            )
            post_slope_residual = p_right - (d_h + e_h - 0.5 * jump)
            result["maximum_post_slope_identity_residual"] = max(
                result["maximum_post_slope_identity_residual"],
                float(np.max(np.abs(post_slope_residual))),
            )
            result["event_count"] += int(q.size)
            add_sign_counts(result["E_h_sign_counts"], e_h)
            add_sign_counts(result["M_T_sign_counts"], reserve)
            result["raw_negative_U_count"] += int(
                np.count_nonzero(turan_defect < 0.0)
            )
            result["U_below_minus_1e_minus_12_count"] += int(
                np.count_nonzero(turan_defect < -1e-12)
            )
            result["nonpositive_phi_2h_count"] += int(
                np.count_nonzero(phi_2h <= 0.0)
            )
            result["nonpositive_radius_squared_count"] += int(
                np.count_nonzero(radius_squared <= 0.0)
            )

            def row_builder(index: int) -> dict[str, Any]:
                return compact_row(
                    index,
                    q,
                    prime,
                    exponent,
                    lambda_value,
                    h,
                    phi,
                    phi_minus,
                    phi_plus,
                    phi_2h,
                    p_bar,
                    d_h,
                    e_h,
                    radius,
                    jump,
                    turan_bound,
                    turan_defect,
                    reserve,
                    actual_post_reserve,
                )

            normalized_u = turan_defect / np.maximum(4.0 * phi_product, 1e-300)
            metric_specs = (
                ("minimum_M_T", reserve, np.argmin, False),
                ("minimum_M_T_over_J", reserve / jump, np.argmin, False),
                ("minimum_actual_post_reserve", actual_post_reserve, np.argmin, False),
                ("minimum_conditional_gap", conditional_gap, np.argmin, False),
                (
                    "minimum_R_minus_turan_bound",
                    radius - turan_bound,
                    np.argmin,
                    False,
                ),
                (
                    "minimum_relative_R_minus_turan_bound",
                    (radius - turan_bound) / radius,
                    np.argmin,
                    False,
                ),
                (
                    "minimum_M_T_over_R",
                    reserve / radius,
                    np.argmin,
                    False,
                ),
                (
                    "maximum_turan_bound_over_R",
                    turan_bound / radius,
                    np.argmax,
                    True,
                ),
                ("minimum_U", turan_defect, np.argmin, False),
                ("minimum_normalized_U", normalized_u, np.argmin, False),
                ("minimum_E_h", e_h, np.argmin, False),
                ("maximum_E_h", e_h, np.argmax, True),
                ("minimum_E_h_over_J", e_h / jump, np.argmin, False),
                ("maximum_E_h_over_J", e_h / jump, np.argmax, True),
                ("maximum_abs_E_h", np.abs(e_h), np.argmax, True),
                ("maximum_abs_E_h_over_J", np.abs(e_h / jump), np.argmax, True),
                (
                    "maximum_abs_E_h_over_eta",
                    np.abs(e_h / (log_x / root_x)),
                    np.argmax,
                    True,
                ),
                ("maximum_abs_E_h_over_R", np.abs(e_h / radius), np.argmax, True),
                (
                    "maximum_abs_E_h_over_turan_bound",
                    np.abs(e_h / turan_bound),
                    np.argmax,
                    True,
                ),
            )
            for name, values, selector, maximum in metric_specs:
                index = int(selector(values))
                update_extreme(
                    result["global_extrema"],
                    name,
                    float(values[index]),
                    row_builder(index),
                    maximum=maximum,
                )
            nonpositive_m_positions = np.flatnonzero(reserve <= 0.0)
            if nonpositive_m_positions.size:
                index = int(nonpositive_m_positions[-1])
                update_extreme(
                    result["global_extrema"],
                    "maximum_q_with_nonpositive_M_T",
                    float(q[index]),
                    row_builder(index),
                    maximum=True,
                )
            negative_u_positions = np.flatnonzero(turan_defect < 0.0)
            if negative_u_positions.size:
                index = int(negative_u_positions[-1])
                update_extreme(
                    result["global_extrema"],
                    "maximum_q_with_negative_U",
                    float(q[index]),
                    row_builder(index),
                    maximum=True,
                )

            for raw_decade in np.unique(decade_keys):
                decade = int(raw_decade)
                positions = np.flatnonzero(decade_keys == raw_decade)
                record = result["decades"].setdefault(
                    decade, new_decade(decade)
                )
                record["event_count"] += int(positions.size)
                add_sign_counts(record["E_h_sign_counts"], e_h[positions])
                add_sign_counts(record["M_T_sign_counts"], reserve[positions])
                record["raw_negative_U_count"] += int(
                    np.count_nonzero(turan_defect[positions] < 0.0)
                )
                record["U_below_minus_1e_minus_12_count"] += int(
                    np.count_nonzero(turan_defect[positions] < -1e-12)
                )
                for name, values, selector, maximum in (
                    ("minimum_M_T", reserve, np.argmin, False),
                    ("minimum_M_T_over_J", reserve / jump, np.argmin, False),
                    ("minimum_U", turan_defect, np.argmin, False),
                    (
                        "minimum_R_minus_turan_bound",
                        radius - turan_bound,
                        np.argmin,
                        False,
                    ),
                    (
                        "minimum_relative_R_minus_turan_bound",
                        (radius - turan_bound) / radius,
                        np.argmin,
                        False,
                    ),
                    (
                        "minimum_M_T_over_R",
                        reserve / radius,
                        np.argmin,
                        False,
                    ),
                    (
                        "maximum_turan_bound_over_R",
                        turan_bound / radius,
                        np.argmax,
                        True,
                    ),
                    ("minimum_normalized_U", normalized_u, np.argmin, False),
                    ("minimum_E_h_over_J", e_h / jump, np.argmin, False),
                    ("maximum_E_h_over_J", e_h / jump, np.argmax, True),
                    ("maximum_abs_E_h", np.abs(e_h), np.argmax, True),
                    (
                        "maximum_abs_E_h_over_J",
                        np.abs(e_h / jump),
                        np.argmax,
                        True,
                    ),
                    (
                        "maximum_abs_E_h_over_R",
                        np.abs(e_h / radius),
                        np.argmax,
                        True,
                    ),
                ):
                    index = int(positions[int(selector(values[positions]))])
                    update_extreme(
                        record["extrema"],
                        name,
                        float(values[index]),
                        row_builder(index),
                        maximum=maximum,
                    )

        processed_chunks += 1
        if args.progress_every and (
            processed_chunks % args.progress_every == 0 or stop == center_count
        ):
            summaries = ", ".join(
                f"lambda={value:g}:minM={lambda_results[lambda_key(value)]['global_extrema']['minimum_M_T']['value']:.6g}"
                for value in lambda_values
            )
            print(
                f"chunk {processed_chunks}: events={stop:,}/{center_count:,}; {summaries}",
                flush=True,
            )

    finalized_results: dict[str, Any] = {}
    for key, result in lambda_results.items():
        decade_rows = [
            result["decades"][decade] for decade in sorted(result["decades"])
        ]
        fits: list[dict[str, Any]] = []
        for metric_name, slope_name in (
            ("maximum_abs_E_h", "slope_log_abs_E_on_log_q"),
            (
                "maximum_abs_E_h_over_R",
                "slope_log_abs_E_over_R_on_log_q",
            ),
        ):
            for cutoff in (10**2, 10**4, 10**6):
                selected = [
                    row["extrema"][metric_name]
                    for row in decade_rows
                    if row["q_right"] >= cutoff
                    and row["extrema"][metric_name]["row"]["q"] >= cutoff
                ]
                if len(selected) < 3:
                    continue
                log_q = np.log(
                    np.asarray(
                        [record["row"]["q"] for record in selected],
                        dtype=np.float64,
                    )
                )
                log_metric = np.log(
                    np.asarray([record["value"] for record in selected])
                )
                slope, intercept = np.polyfit(log_q, log_metric, 1)
                predicted = intercept + slope * log_q
                residual = float(np.sum((log_metric - predicted) ** 2))
                total = float(
                    np.sum((log_metric - float(np.mean(log_metric))) ** 2)
                )
                fits.append(
                    {
                        "metric": metric_name,
                        "cutoff": cutoff,
                        "decade_count": len(selected),
                        slope_name: float(slope),
                        "intercept": float(intercept),
                        "R_squared": 1.0 - residual / total if total else None,
                    }
                )
        result["decade_diagnostics"] = decade_rows
        del result["decades"]
        result["effective_scaling_fits"] = fits
        result["checks"] = {
            "all_phi_2h_positive": result["nonpositive_phi_2h_count"] == 0,
            "all_radius_squared_positive": (
                result["nonpositive_radius_squared_count"] == 0
            ),
            "algebra_identity_residual_below_1e_minus_10": (
                result["maximum_algebra_identity_residual"] < 1e-10
            ),
            "post_slope_identity_residual_below_1e_minus_12": (
                result["maximum_post_slope_identity_residual"] < 1e-12
            ),
        }
        result["checks"]["all_internal_checks_pass"] = all(
            result["checks"].values()
        )
        finalized_results[key] = result

    support_prime_count = int(np.count_nonzero(event_exponent == 1))
    support_higher_count = int(event_q.size - support_prime_count)
    center_prime_count = int(np.count_nonzero(event_exponent[:center_count] == 1))
    center_higher_count = int(center_count - center_prime_count)
    report = {
        "status": (
            "Exploratory binary64 prime-side scan. U is an RH-conditional "
            "Turan diagnostic; finite nonnegativity is not a proof of RH."
        ),
        "configuration": {
            "center_limit": args.limit,
            "support_limit": support_limit,
            "segment_span": args.segment_span,
            "chunk_size": args.chunk_size,
            "lambda_values": list(lambda_values),
        },
        "definitions": {
            "h": "lambda*log(q)/sqrt(q)",
            "D_h": "(Phi(t+h)-Phi(t-h))/(2h)",
            "p_bar": "(Phi_prime(t-)+Phi_prime(t+))/2",
            "E_h": "p_bar-D_h",
            "U": "4*Phi(t)*Phi(2h)-(Phi(t+h)-Phi(t-h))^2",
            "M_T": "R-sqrt(Phi(t)*Phi(2h))/h+E_h-J/2",
            "actual_post_reserve": "R+p_bar-J/2=R+Phi_prime(t+)",
            "conditional_implication": (
                "If U>=0 then actual_post_reserve>=M_T; hence M_T>0 "
                "is a sufficient RH-conditional event test."
            ),
        },
        "counts": {
            "center_prime_events": center_prime_count,
            "center_higher_prime_power_events": center_higher_count,
            "center_all_prime_power_events": center_count,
            "support_prime_events": support_prime_count,
            "support_higher_prime_power_events": support_higher_count,
            "support_all_prime_power_events": int(event_q.size),
        },
        "lambda_results": finalized_results,
        "checks": {
            "center_event_count_decomposes": (
                center_count == center_prime_count + center_higher_count
            ),
            "support_event_count_decomposes": (
                event_q.size == support_prime_count + support_higher_count
            ),
            "all_lambda_internal_checks_pass": all(
                result["checks"]["all_internal_checks_pass"]
                for result in finalized_results.values()
            ),
        },
        "runtime": {
            "seconds": time.perf_counter() - started,
            "python": platform.python_version(),
            "numpy": np.__version__,
            "platform": platform.platform(),
            "script_sha256": sha256(SCRIPT_PATH),
            "helper_sha256": sha256(HELPER_PATH),
        },
    }
    report["checks"]["all_pass"] = all(report["checks"].values())
    return report


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=10**7)
    parser.add_argument("--segment-span", type=int, default=10_000_000)
    parser.add_argument("--chunk-size", type=int, default=250_000)
    parser.add_argument(
        "--lambdas", type=float, nargs="+", default=list(DEFAULT_LAMBDAS)
    )
    parser.add_argument("--progress-every", type=int, default=4)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def main() -> None:
    """Run the scan and write its JSON report."""

    args = parse_args()
    if args.limit < 2:
        raise SystemExit("--limit must be at least 2")
    if args.segment_span < 1 or args.chunk_size < 1:
        raise SystemExit("--segment-span and --chunk-size must be positive")
    if not args.lambdas or any(
        value <= 0.0 or value >= math.sqrt(2.0) for value in args.lambdas
    ):
        raise SystemExit(
            "all lambda values must satisfy 0<lambda<sqrt(2) so every window has 0<h<t"
        )
    report = scan(args)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {args.report}")
    print(
        json.dumps(
            {
                "center_limit": args.limit,
                "event_count": report["counts"]["center_all_prime_power_events"],
                "lambda_summary": {
                    key: {
                        "minimum_M_T": value["global_extrema"]["minimum_M_T"],
                        "minimum_U": value["global_extrema"]["minimum_U"],
                    }
                    for key, value in report["lambda_results"].items()
                },
                "all_checks_pass": report["checks"]["all_pass"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
