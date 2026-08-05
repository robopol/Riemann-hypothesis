#!/usr/bin/env python3
"""Explore exact prime-side forms of the integrated Robin residual.

This script is an exploratory binary64 computation, not an interval
certificate and not a proof of the Riemann hypothesis.  It performs five
related checks:

1. reconstruct R_infinity from a finite prime-power sum and compare it with
   the independent E/C_pp ledger and finite block integrals;
2. compare the direct D_* prime sum with its exact lead-window form and scan
   the modified Chebyshev signal Psi_*(x)-x at prime supports;
3. evaluate the prime-side Weil function B(log x) and compare it with an
   Odlyzko on-line-zero partial sum;
4. form several sampled Gram/Toeplitz matrices B(t_j-t_k);
5. scan the necessary positive-definite bound |B(t)| <= B(0).

All code and report labels explicitly distinguish exact symbolic identities
from floating-point numerical evidence.
"""

from __future__ import annotations

import argparse
import bisect
import hashlib
import json
import math
import platform
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np


ROOT = Path(__file__).resolve().parent
DEFAULT_ZEROS = ROOT / "exploratory_integral_zeros1.txt"
DEFAULT_REPORT = ROOT / "exploratory_integral_prime_side_report.json"

EULER_GAMMA = 0.577215664901532860606512090082402431
LOG_TWO_PI = math.log(2.0 * math.pi)
WEIL_C = 2.0 + EULER_GAMMA - math.log(4.0 * math.pi)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def odd_sieve(limit: int) -> list[int]:
    """Return all primes through limit with a deterministic odd-only sieve."""

    if limit < 2:
        return []
    flags = bytearray(b"\x01") * (limit // 2 + 1)
    flags[0] = 0
    for prime in range(3, math.isqrt(limit) + 1, 2):
        if flags[prime // 2]:
            start = prime * prime // 2
            count = (len(flags) - start - 1) // prime + 1
            flags[start::prime] = b"\x00" * count
    return [2] + [value for value in range(3, limit + 1, 2) if flags[value // 2]]


@dataclass(frozen=True)
class PrimeTables:
    """Prefix tables for primes and prime-power events."""

    limit: int
    primes: np.ndarray
    prime_log: np.ndarray
    theta_prefix: np.ndarray
    log_beta_prefix: np.ndarray
    event_value: np.ndarray
    event_log_prime: np.ndarray
    psi_prefix: np.ndarray
    s_log_prefix: np.ndarray
    s_lambda_prefix: np.ndarray

    def prime_count(self, x: float, side: Literal["left", "right"] = "right") -> int:
        return int(np.searchsorted(self.primes, x, side=side))

    def event_count(self, x: float, side: Literal["left", "right"] = "right") -> int:
        return int(np.searchsorted(self.event_value, x, side=side))

    @staticmethod
    def _prefix_value(prefix: np.ndarray, count: int) -> float:
        return float(prefix[count - 1]) if count else 0.0

    def query(
        self, x: float, side: Literal["left", "right"] = "right"
    ) -> dict[str, float | int]:
        prime_count = self.prime_count(x, side)
        event_count = self.event_count(x, side)
        return {
            "prime_count": prime_count,
            "event_count": event_count,
            "theta": self._prefix_value(self.theta_prefix, prime_count),
            "log_beta": self._prefix_value(self.log_beta_prefix, prime_count),
            "psi": self._prefix_value(self.psi_prefix, event_count),
            "s_log": self._prefix_value(self.s_log_prefix, event_count),
            "s_lambda": self._prefix_value(self.s_lambda_prefix, event_count),
        }


def build_prime_tables(limit: int) -> PrimeTables:
    """Build prime and prime-power prefix sums through limit."""

    primes_list = odd_sieve(limit)
    primes = np.asarray(primes_list, dtype=np.int64)
    prime_log = np.log(primes.astype(np.float64))
    theta_prefix = np.cumsum(prime_log, dtype=np.float64)
    log_beta_prefix = np.cumsum(
        -np.log1p(-1.0 / primes.astype(np.float64)), dtype=np.float64
    )

    events: list[tuple[int, float, float, float]] = []
    for prime, log_prime in zip(primes_list, prime_log, strict=True):
        value = prime
        exponent = 1
        while value <= limit:
            events.append(
                (
                    value,
                    float(log_prime),
                    1.0 / (exponent * value),
                    float(log_prime) / value,
                )
            )
            if value > limit // prime:
                break
            value *= prime
            exponent += 1
    events.sort(key=lambda row: row[0])

    event_value = np.fromiter((row[0] for row in events), dtype=np.int64)
    event_log_prime = np.fromiter((row[1] for row in events), dtype=np.float64)
    event_s_log = np.fromiter((row[2] for row in events), dtype=np.float64)
    event_s_lambda = np.fromiter((row[3] for row in events), dtype=np.float64)
    return PrimeTables(
        limit=limit,
        primes=primes,
        prime_log=prime_log,
        theta_prefix=theta_prefix,
        log_beta_prefix=log_beta_prefix,
        event_value=event_value,
        event_log_prime=event_log_prime,
        psi_prefix=np.cumsum(event_log_prime, dtype=np.float64),
        s_log_prefix=np.cumsum(event_s_log, dtype=np.float64),
        s_lambda_prefix=np.cumsum(event_s_lambda, dtype=np.float64),
    )


def reciprocal_kernel(x: float) -> float:
    return 1.0 / (x * math.log(x))


def rinf_finite_prime_side(
    tables: PrimeTables, x: float, side: Literal["left", "right"] = "right"
) -> float:
    """Evaluate the exact finite prime-side identity for R_infinity."""

    state = tables.query(x, side)
    log_x = math.log(x)
    return (
        (float(state["psi"]) - x) / (x * log_x)
        + EULER_GAMMA
        + math.log(log_x)
        - float(state["s_log"])
    )


def rinf_ledger_form(
    tables: PrimeTables, x: float, side: Literal["left", "right"] = "right"
) -> float:
    """Evaluate R_infinity as C_pp + endpoint - E."""

    state = tables.query(x, side)
    log_x = math.log(x)
    scale = x * log_x
    cross_tail = float(state["log_beta"]) - float(state["s_log"])
    c_pp = (
        (float(state["psi"]) - float(state["theta"])) / scale + cross_tail
    )
    mertens_error = (
        float(state["log_beta"]) - EULER_GAMMA - math.log(log_x)
    )
    endpoint = (float(state["theta"]) - x) / scale
    return c_pp + endpoint - mertens_error


def rinf_block_integral(tables: PrimeTables, x: float, y: float) -> float:
    """Integrate (psi(t)-t) omega(t) exactly over a finite block."""

    if not 1.0 < x < y <= tables.limit:
        raise ValueError("block endpoints must satisfy 1 < x < y <= table limit")
    x_state = tables.query(x, "right")
    start = tables.event_count(x, "right")
    stop = tables.event_count(y, "right")
    f_x = reciprocal_kernel(x)
    f_y = reciprocal_kernel(y)
    psi_integral = float(x_state["psi"]) * (f_x - f_y)
    if stop > start:
        values = tables.event_value[start:stop].astype(np.float64)
        logs = tables.event_log_prime[start:stop]
        psi_integral += float(np.sum(logs * (1.0 / (values * np.log(values)) - f_y)))
    main_integral = (
        math.log(math.log(y) / math.log(x))
        + 1.0 / math.log(x)
        - 1.0 / math.log(y)
    )
    return psi_integral - main_integral


def rinf_checks(tables: PrimeTables) -> dict[str, Any]:
    points = [
        3_299.0,
        4_391.0,
        10_007.0,
        87_833.0,
        319_439.0,
        739_951.0,
        1_412_351.0,
        3_329_267.0,
    ]
    points = [x for x in points if x <= tables.limit]
    rows: list[dict[str, float]] = []
    for x in points:
        finite = rinf_finite_prime_side(tables, x)
        ledger = rinf_ledger_form(tables, x)
        rows.append(
            {
                "x": x,
                "finite_prime_side": finite,
                "ledger_cpp_endpoint_e": ledger,
                "difference": finite - ledger,
                "normalized_rinf": math.sqrt(x) * math.log(x) * finite,
            }
        )

    block_pairs = [
        (3_299.25, 10_007.25),
        (10_007.25, 87_833.25),
        (87_833.25, 739_951.25),
        (739_951.25, 3_329_267.25),
    ]
    block_rows: list[dict[str, float]] = []
    for x, y in block_pairs:
        if y > tables.limit:
            continue
        direct = rinf_block_integral(tables, x, y)
        tail_difference = rinf_finite_prime_side(tables, x) - rinf_finite_prime_side(
            tables, y
        )
        block_rows.append(
            {
                "x": x,
                "y": y,
                "direct_finite_block_integral": direct,
                "rinf_x_minus_rinf_y": tail_difference,
                "difference": direct - tail_difference,
            }
        )

    endpoint_rows: list[dict[str, float]] = []
    for value in (2, 4, 8, 9, 27, 64, 125, 729):
        if value > tables.limit:
            continue
        left = rinf_finite_prime_side(tables, float(value), "left")
        right = rinf_finite_prime_side(tables, float(value), "right")
        endpoint_rows.append(
            {
                "prime_power": float(value),
                "left_convention": left,
                "right_convention": right,
                "difference": right - left,
            }
        )

    return {
        "identity": (
            "R_inf(x)=(psi(x)-x)/(x log x)+gamma_E+log log x"
            "-sum_{p^m<=x}1/(m p^m)"
        ),
        "selected_points": rows,
        "finite_block_checks": block_rows,
        "prime_power_endpoint_continuity": endpoint_rows,
        "max_abs_ledger_difference": max(abs(row["difference"]) for row in rows),
        "max_abs_block_difference": max(
            abs(row["difference"]) for row in block_rows
        ),
        "max_abs_endpoint_difference": max(
            abs(row["difference"]) for row in endpoint_rows
        ),
    }


def optimizer_delta(prime: int, layer: int) -> float:
    """Return the positive A-drop when exponent layer-1 rises to layer."""

    inverse_power = float(prime) ** (-layer)
    return math.log1p(-inverse_power / prime) - math.log1p(-inverse_power)


def prime_tail_q(prime: int, k: int, tolerance: float = 1e-24) -> float:
    value = prime ** (k + 1)
    exponent = k + 1
    total = 0.0
    while True:
        term = 1.0 / (exponent * value)
        total += term
        if term < tolerance:
            return total
        value *= prime
        exponent += 1


def dstar_direct_and_windows(
    x: int, primes: np.ndarray, tolerance: float = 1e-22
) -> dict[str, float | int]:
    """Compare direct and lead-window forms of D_* at one support."""

    log_x = math.log(x)
    f_x = 1.0 / (x * log_x)
    direct = 0.0
    windows = 0.0
    maximum_optimizer_exponent = 1
    for prime_value in primes:
        prime = int(prime_value)
        if prime > x:
            break
        log_prime = math.log(prime)
        c = log_prime * f_x

        k = 1
        natural_power = prime
        while natural_power <= x // prime:
            natural_power *= prime
            k += 1

        optimizer_exponent = 1
        inverse_power = 1.0 / (prime * prime)
        delta = math.log1p(-inverse_power / prime) - math.log1p(-inverse_power)
        while delta > c:
            optimizer_exponent += 1
            inverse_power /= prime
            delta = math.log1p(-inverse_power / prime) - math.log1p(-inverse_power)
        maximum_optimizer_exponent = max(
            maximum_optimizer_exponent, optimizer_exponent
        )
        d_value = -math.log1p(-inverse_power)
        direct += (
            (optimizer_exponent - k) * c
            + d_value
            - prime_tail_q(prime, k)
        )

        layer = 2
        layer_power = prime * prime
        inverse_layer_power = 1.0 / layer_power
        while True:
            layer_delta = math.log1p(
                -inverse_layer_power / prime
            ) - math.log1p(-inverse_layer_power)
            natural_mass = 1.0 / (layer * layer_power)
            if layer_power > x:
                windows += min(c, layer_delta) - natural_mass
            if max(layer_delta, natural_mass) < tolerance:
                break
            layer += 1
            layer_power *= prime
            inverse_layer_power /= prime

    scale = math.sqrt(x) * log_x
    return {
        "x": x,
        "direct_dstar": direct,
        "lead_window_dstar": windows,
        "difference": direct - windows,
        "scaled_direct_dstar": scale * direct,
        "scaled_lead_window_dstar": scale * windows,
        "maximum_optimizer_exponent": maximum_optimizer_exponent,
    }


def psi_star_scan(
    tables: PrimeTables, scan_min: int, scan_limit: int
) -> dict[str, Any]:
    """Scan Psi_*(x)-x at every prime support in the requested range."""

    stop = int(np.searchsorted(tables.primes, scan_limit, side="right"))
    primes = tables.primes[:stop]
    logs = tables.prime_log[:stop]
    controls = 1.0 / (primes.astype(np.float64) * logs)
    optimizer_add = np.zeros(stop, dtype=np.float64)
    event_count = 0
    final_control = float(controls[-1])

    candidate_cutoff = math.isqrt(2 * scan_limit) + 2
    for prime_value, log_prime in zip(primes, logs, strict=True):
        prime = int(prime_value)
        if prime > candidate_cutoff:
            break
        layer = 2
        inverse_power = 1.0 / (prime * prime)
        while True:
            delta = math.log1p(-inverse_power / prime) - math.log1p(-inverse_power)
            threshold = delta / float(log_prime)
            if threshold <= final_control:
                break
            activation = int(np.searchsorted(-controls, -threshold, side="right"))
            if activation < stop:
                optimizer_add[activation] += float(log_prime)
                event_count += 1
            layer += 1
            inverse_power /= prime

    power_events: list[tuple[int, float]] = []
    for prime_value, log_prime in zip(primes, logs, strict=True):
        prime = int(prime_value)
        if prime * prime > scan_limit:
            break
        value = prime * prime
        while value <= scan_limit:
            power_events.append((value, float(log_prime)))
            if value > scan_limit // prime:
                break
            value *= prime
    power_events.sort(key=lambda row: row[0])

    theta = 0.0
    optimizer_h = 0.0
    psi_minus_theta = 0.0
    power_index = 0
    support_count = 0
    nonpositive_count = 0
    minimum: dict[str, float] | None = None
    maximum: dict[str, float] | None = None
    first_nonpositive: dict[str, float] | None = None
    last_nonpositive: dict[str, float] | None = None
    endpoint: dict[str, float] | None = None

    for index, (prime_value, log_prime) in enumerate(zip(primes, logs, strict=True)):
        prime = int(prime_value)
        theta += float(log_prime)
        optimizer_h += float(optimizer_add[index])
        while power_index < len(power_events) and power_events[power_index][0] <= prime:
            psi_minus_theta += power_events[power_index][1]
            power_index += 1
        if prime < scan_min:
            continue
        support_count += 1
        value = theta + optimizer_h - prime
        scaled = value / math.sqrt(prime)
        row = {
            "support": float(prime),
            "psi_star_minus_x": value,
            "scaled_by_sqrt_x": scaled,
            "psi_minus_x": theta + psi_minus_theta - prime,
            "lead_window_mass_K": optimizer_h - psi_minus_theta,
        }
        if minimum is None or scaled < minimum["scaled_by_sqrt_x"]:
            minimum = row.copy()
        if maximum is None or scaled > maximum["scaled_by_sqrt_x"]:
            maximum = row.copy()
        if value <= 0.0:
            nonpositive_count += 1
            if first_nonpositive is None:
                first_nonpositive = row.copy()
            last_nonpositive = row.copy()
        endpoint = row.copy()

    return {
        "status": (
            "Binary64 threshold scan at every prime support; exploratory only, "
            "not an interval certificate."
        ),
        "definition": "Psi_*(x)=log n_*(x)=psi(x)+K_*(x)",
        "scan_min": scan_min,
        "scan_limit": scan_limit,
        "support_count": support_count,
        "optimizer_activation_event_count": event_count,
        "prime_power_event_count": len(power_events),
        "nonpositive_count": nonpositive_count,
        "minimum": minimum,
        "maximum": maximum,
        "first_nonpositive": first_nonpositive,
        "last_nonpositive": last_nonpositive,
        "endpoint": endpoint,
    }


def lead_window_checks(
    tables: PrimeTables, scan_min: int, scan_limit: int
) -> dict[str, Any]:
    points = [3_299, 10_007, 87_833, 739_951]
    points = [x for x in points if x <= scan_limit]
    rows = [dstar_direct_and_windows(x, tables.primes) for x in points]
    return {
        "identity": (
            "D_*(x)=integral_x^infinity K_*(t) omega(t) dt"
            "-sum_{p>x}(1/p-log(1+1/p))"
        ),
        "combined_tail_identity": (
            "R_inf(x)+D_*(x)=integral_x^infinity (Psi_*(t)-t) omega(t) dt"
            "-sum_{p>x}(1/p-log(1+1/p))"
        ),
        "selected_dstar_checks": rows,
        "max_abs_direct_minus_windows": max(abs(float(row["difference"])) for row in rows),
        "psi_star_scan": psi_star_scan(tables, scan_min, scan_limit),
    }


def trivial_f(y: float) -> float:
    """Return F(y)=2 atanh(y)-2y+y log(1-y^2) stably."""

    if not 0.0 < y < 1.0:
        raise ValueError("trivial-factor argument must lie in (0,1)")
    if y < 0.25:
        # F(y)=-sum_{k>=1} y^(2k+1)/(k(2k+1)).
        power = y * y * y
        total = 0.0
        layer = 1
        while True:
            term = -power / (layer * (2 * layer + 1))
            total += term
            if abs(term) < 1e-24:
                return total
            power *= y * y
            layer += 1
    return (
        (1.0 + y) * math.log1p(y)
        + (y - 1.0) * math.log1p(-y)
        - 2.0 * y
    )


def prime_side_weil_b(
    tables: PrimeTables,
    t: float,
    side: Literal["left", "right"] = "right",
) -> float:
    """Evaluate the prime-side Weil function B(t), using x=exp(|t|)."""

    t = abs(t)
    if t == 0.0:
        return WEIL_C
    x = math.exp(t)
    if x > tables.limit:
        raise ValueError("Weil evaluation exceeds the prime table limit")
    state = tables.query(x, side)
    y = 1.0 / x
    bracket = (
        (float(state["s_lambda"]) - t + EULER_GAMMA)
        - (float(state["psi"]) / x - 1.0)
        - LOG_TWO_PI / x
        - 0.5 * trivial_f(y)
    )
    return math.sqrt(x) * bracket


def load_zeros(path: Path, count: int) -> np.ndarray:
    values = np.loadtxt(path, dtype=np.float64, max_rows=count)
    if values.ndim != 1 or values.size == 0:
        raise ValueError(f"no zero ordinates found in {path}")
    if not np.all(np.diff(values) > 0.0):
        raise ValueError("zero ordinates are not strictly increasing")
    return values


def zero_partial_weil_b(gamma: np.ndarray, t_values: np.ndarray) -> np.ndarray:
    """Evaluate 2 sum cos(gamma*t)/(gamma^2+1/4) in bounded chunks."""

    weights = 2.0 / (gamma * gamma + 0.25)
    output = np.empty(t_values.size, dtype=np.float64)
    chunk_size = 32
    for start in range(0, t_values.size, chunk_size):
        stop = min(start + chunk_size, t_values.size)
        phase = np.outer(t_values[start:stop], gamma)
        output[start:stop] = np.cos(phase) @ weights
    return output


def weil_comparison(
    tables: PrimeTables,
    gamma: np.ndarray,
    t_max: float,
    grid_size: int,
) -> dict[str, Any]:
    grid = np.linspace(0.0, t_max, grid_size)
    prime_values = np.asarray(
        [prime_side_weil_b(tables, float(t)) for t in grid], dtype=np.float64
    )
    zero_values = zero_partial_weil_b(gamma, grid)
    differences = prime_values - zero_values
    partial_mass = float(np.sum(2.0 / (gamma * gamma + 0.25)))

    selected_t = [0.0, 0.05, 0.1, 0.5, 0.75, 1.0, 2.0, 5.0, 10.0, t_max]
    selected_t = sorted(set(t for t in selected_t if t <= t_max))
    selected_zero = zero_partial_weil_b(gamma, np.asarray(selected_t, dtype=np.float64))
    selected_rows = []
    for t, zero_value in zip(selected_t, selected_zero, strict=True):
        prime_value = prime_side_weil_b(tables, t)
        selected_rows.append(
            {
                "t": t,
                "x": math.exp(t),
                "prime_side_B": prime_value,
                "odlyzko_partial_B": float(zero_value),
                "prime_minus_partial": prime_value - float(zero_value),
            }
        )

    endpoint_rows = []
    for value in (2, 4, 8, 9, 27, 64, 125, 729):
        if value > tables.limit:
            continue
        t = math.log(value)
        left = prime_side_weil_b(tables, t, "left")
        right = prime_side_weil_b(tables, t, "right")
        endpoint_rows.append(
            {
                "prime_power": value,
                "t": t,
                "left_convention": left,
                "right_convention": right,
                "difference": right - left,
            }
        )

    return {
        "identity": (
            "B(log x)=sqrt(x)[(S_Lambda(x)-log x+gamma_E)"
            "-(psi(x)/x-1)-log(2pi)/x-F(1/x)/2]"
        ),
        "definitions": {
            "S_Lambda": "sum_{n<=x} Lambda(n)/n",
            "F": "F(y)=2 atanh(y)-2y+y log(1-y^2)",
            "zero_side": "B(t)=sum_rho exp((rho-1/2)t)/(rho(1-rho))",
            "on_RH": "B(t)=2 sum_{gamma>0} cos(gamma t)/(gamma^2+1/4)",
        },
        "zero_count": int(gamma.size),
        "partial_zero_mass_at_t0": partial_mass,
        "exact_B0": WEIL_C,
        "unaccounted_B0_mass": WEIL_C - partial_mass,
        "grid": {
            "t_min": 0.0,
            "t_max": t_max,
            "size": grid_size,
            "max_abs_prime_minus_partial": float(np.max(np.abs(differences))),
            "max_abs_prime_minus_partial_excluding_t0": float(
                np.max(np.abs(differences[1:]))
            ),
            "rms_prime_minus_partial": float(np.sqrt(np.mean(differences**2))),
        },
        "selected_points": selected_rows,
        "prime_power_endpoint_continuity": endpoint_rows,
        "max_abs_endpoint_difference": max(
            abs(row["difference"]) for row in endpoint_rows
        ),
    }


def gram_matrix_report(tables: PrimeTables) -> dict[str, Any]:
    specifications: list[tuple[str, list[float]]] = [
        ("toeplitz_n4_step_0p75", [0.75 * index for index in range(4)]),
        ("toeplitz_n8_step_1", [float(index) for index in range(8)]),
        ("toeplitz_n12_step_1", [float(index) for index in range(12)]),
        ("irregular_n7", [0.0, 0.55, 1.4, 2.8, 4.7, 7.5, 10.5]),
    ]
    rows: list[dict[str, Any]] = []
    for name, points in specifications:
        matrix = np.asarray(
            [
                [prime_side_weil_b(tables, left - right) for right in points]
                for left in points
            ],
            dtype=np.float64,
        )
        eigenvalues = np.linalg.eigvalsh(matrix)
        rows.append(
            {
                "name": name,
                "points": points,
                "matrix": matrix.tolist(),
                "eigenvalues": eigenvalues.tolist(),
                "minimum_eigenvalue": float(eigenvalues[0]),
                "maximum_eigenvalue": float(eigenvalues[-1]),
                "symmetry_max_abs_error": float(np.max(np.abs(matrix - matrix.T))),
                "classification": (
                    "sampled_positive_definite"
                    if eigenvalues[0] > 0.0
                    else "sampled_not_positive_definite"
                ),
            }
        )
    return {
        "status": (
            "Floating sampled matrices only. Positive sampled eigenvalues are "
            "not a proof of positive definiteness or RH."
        ),
        "rows": rows,
        "smallest_minimum_eigenvalue": min(
            row["minimum_eigenvalue"] for row in rows
        ),
    }


def weil_bound_scan(
    tables: PrimeTables, t_max: float, grid_size: int
) -> dict[str, Any]:
    grid = np.linspace(0.0, t_max, grid_size)
    values = np.asarray(
        [prime_side_weil_b(tables, float(t)) for t in grid], dtype=np.float64
    )
    absolute = np.abs(values)
    all_index = int(np.argmax(absolute))
    nonzero_index = int(np.argmax(absolute[1:])) + 1
    minimum_index = int(np.argmin(values))
    maximum_index = int(np.argmax(values))
    return {
        "status": (
            "Floating finite-grid test of the necessary positive-definite bound; "
            "not a uniform theorem."
        ),
        "bound_C_equals_B0": WEIL_C,
        "t_min": 0.0,
        "t_max": t_max,
        "grid_size": grid_size,
        "maximum_absolute_including_t0": {
            "t": float(grid[all_index]),
            "value": float(values[all_index]),
            "absolute_value": float(absolute[all_index]),
            "C_minus_absolute_value": WEIL_C - float(absolute[all_index]),
        },
        "maximum_absolute_excluding_t0": {
            "t": float(grid[nonzero_index]),
            "value": float(values[nonzero_index]),
            "absolute_value": float(absolute[nonzero_index]),
            "C_minus_absolute_value": WEIL_C - float(absolute[nonzero_index]),
        },
        "minimum_B": {
            "t": float(grid[minimum_index]),
            "value": float(values[minimum_index]),
        },
        "maximum_B": {
            "t": float(grid[maximum_index]),
            "value": float(values[maximum_index]),
        },
        "sampled_bound_holds_with_1e_minus_10_tolerance": bool(
            np.max(absolute) <= WEIL_C + 1e-10
        ),
    }


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    started = time.perf_counter()
    required_limit = max(
        args.prime_limit,
        args.scan_limit,
        math.ceil(math.exp(args.weil_t_max)),
    )
    tables = build_prime_tables(required_limit)
    zeros = load_zeros(args.zeros.resolve(), args.zero_count)

    rinf = rinf_checks(tables)
    lead = lead_window_checks(tables, args.scan_min, args.scan_limit)
    weil = weil_comparison(
        tables, zeros, args.weil_t_max, args.weil_grid_size
    )
    gram = gram_matrix_report(tables)
    bound = weil_bound_scan(tables, args.weil_t_max, args.bound_grid_size)

    checks = {
        "rinf_ledger_cross_check_below_1e_minus_10": (
            rinf["max_abs_ledger_difference"] < 1e-10
        ),
        "rinf_block_cross_check_below_1e_minus_10": (
            rinf["max_abs_block_difference"] < 1e-10
        ),
        "rinf_endpoint_continuity_below_1e_minus_12": (
            rinf["max_abs_endpoint_difference"] < 1e-12
        ),
        "lead_window_identity_below_1e_minus_12": (
            lead["max_abs_direct_minus_windows"] < 1e-12
        ),
        "weil_endpoint_continuity_below_1e_minus_12": (
            weil["max_abs_endpoint_difference"] < 1e-12
        ),
        "prime_vs_100k_zero_partial_below_1e_minus_4": (
            weil["grid"]["max_abs_prime_minus_partial"] < 1e-4
        ),
        "sampled_gram_matrices_have_positive_minimum_eigenvalues": (
            gram["smallest_minimum_eigenvalue"] > 0.0
        ),
        "sampled_abs_B_bound_holds": bound[
            "sampled_bound_holds_with_1e_minus_10_tolerance"
        ],
    }
    checks["all_pass"] = all(checks.values())

    return {
        "status": (
            "Exploratory binary64 computation. Symbolic identities are exact, but "
            "all scans, zero comparisons, Gram eigenvalues, and sampled bounds are "
            "floating-point evidence only; this is not an interval certificate or RH proof."
        ),
        "input": {
            "zero_file": str(args.zeros.resolve()),
            "zero_file_sha256": sha256(args.zeros.resolve()),
            "zero_count": int(zeros.size),
            "first_zero_ordinate": float(zeros[0]),
            "last_zero_ordinate": float(zeros[-1]),
            "prime_table_limit": required_limit,
            "prime_count": int(tables.primes.size),
            "prime_power_event_count": int(tables.event_value.size),
            "scan_min": args.scan_min,
            "scan_limit": args.scan_limit,
            "weil_t_max": args.weil_t_max,
            "weil_grid_size": args.weil_grid_size,
            "bound_grid_size": args.bound_grid_size,
        },
        "formula_sign_audit": {
            "prime_side_B": (
                "+sqrt(x)[(S_Lambda-log x+gamma_E)-(psi/x-1)"
                "-log(2pi)/x-F(1/x)/2]"
            ),
            "zero_side_B": (
                "+sum_rho x^(rho-1/2)/(rho(1-rho)); on RH this is "
                "+2 sum cos(gamma log x)/(gamma^2+1/4)"
            ),
            "endpoint_note": (
                "At x=p^m the jumps Lambda(x)/x in S_Lambda and "
                "Lambda(x)/x in psi(x)/x cancel exactly."
            ),
        },
        "finite_rinf_prime_side": rinf,
        "lead_window_dstar": lead,
        "weil_prime_vs_zeros": weil,
        "gram_toeplitz_matrices": gram,
        "sampled_abs_B_bound": bound,
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
    parser.add_argument("--zeros", type=Path, default=DEFAULT_ZEROS)
    parser.add_argument("--zero-count", type=int, default=100_000)
    parser.add_argument("--prime-limit", type=int, default=5_000_000)
    parser.add_argument("--scan-min", type=int, default=3_299)
    parser.add_argument("--scan-limit", type=int, default=5_000_000)
    parser.add_argument("--weil-t-max", type=float, default=12.0)
    parser.add_argument("--weil-grid-size", type=int, default=241)
    parser.add_argument("--bound-grid-size", type=int, default=2_401)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    if args.scan_min < 2 or args.scan_min > args.scan_limit:
        raise ValueError("scan bounds must satisfy 2 <= scan_min <= scan_limit")
    if args.weil_t_max < 11.0:
        raise ValueError("weil-t-max must be at least 11 for the default Gram grids")
    report = build_report(args)
    if not report["checks"]["all_pass"]:
        raise AssertionError(f"exploratory consistency checks failed: {report['checks']}")
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"wrote {args.report}")
    print(f"all checks pass: {report['checks']['all_pass']}")
    print(
        "lead-window max difference: "
        f"{report['lead_window_dstar']['max_abs_direct_minus_windows']:.3e}"
    )
    print(
        "Weil prime/zero max grid difference: "
        f"{report['weil_prime_vs_zeros']['grid']['max_abs_prime_minus_partial']:.3e}"
    )
    print(
        "smallest sampled Gram eigenvalue: "
        f"{report['gram_toeplitz_matrices']['smallest_minimum_eigenvalue']:.12g}"
    )


if __name__ == "__main__":
    main()
