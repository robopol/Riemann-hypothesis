#!/usr/bin/env python3
"""Scan the eventwise first-crossing reserve inequality far beyond 56 million.

This is an exploratory binary64 computation. It is not an interval
certificate and does not prove the eventwise inequality or RH. The scan uses
an odd-only segmented sieve, merges prime events with all higher prime powers,
and evaluates the exact finite prime-side state immediately before each event.

The main target is

    m(log(q)-) > Lambda(q) / sqrt(q),

where m = R + Phi', Phi = C - B, and
R^2 = Phi^2/4 + 2*k*Phi. In addition to the raw ratio, the script measures
q^(1/4)-normalized lower envelopes and the exact event-to-event reserve
balance m_i^+ = m_{i-1}^+ + smooth_gain_i - J_i.
"""

from __future__ import annotations

import argparse
import hashlib
import heapq
import json
import math
import platform
import time
from pathlib import Path
from typing import Any

import numpy as np


EULER_GAMMA = 0.577215664901532860606512090082402431
LOG_TWO_PI = math.log(2.0 * math.pi)
WEIL_C = 2.0 + EULER_GAMMA - math.log(4.0 * math.pi)
DEFAULT_REPORT = Path(__file__).with_name(
    "exploratory_eventwise_reserve_scan_report.json"
)


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def simple_primes(limit: int) -> list[int]:
    """Return all primes at most limit with a compact deterministic sieve."""

    if limit < 2:
        return []
    flags = bytearray(b"\x01") * (limit + 1)
    flags[0:2] = b"\x00\x00"
    for prime in range(2, math.isqrt(limit) + 1):
        if flags[prime]:
            start = prime * prime
            flags[start : limit + 1 : prime] = b"\x00" * (
                ((limit - start) // prime) + 1
            )
    return [value for value, is_prime in enumerate(flags) if is_prime]


def prime_array_in_segment(
    low: int, high: int, base_primes: list[int]
) -> np.ndarray:
    """Return a sorted int64 array of primes in the inclusive segment."""

    pieces: list[np.ndarray] = []
    if low <= 2 <= high:
        pieces.append(np.asarray([2], dtype=np.int64))

    odd_low = max(3, low)
    if odd_low % 2 == 0:
        odd_low += 1
    if odd_low > high:
        return pieces[0] if pieces else np.empty(0, dtype=np.int64)

    size = ((high - odd_low) // 2) + 1
    flags = np.ones(size, dtype=np.bool_)
    for prime in base_primes[1:]:
        square = prime * prime
        if square > high:
            break
        start = max(square, ((odd_low + prime - 1) // prime) * prime)
        if start % 2 == 0:
            start += prime
        if start > high:
            continue
        flags[(start - odd_low) // 2 :: prime] = False
    odd_primes = odd_low + 2 * np.flatnonzero(flags).astype(np.int64)
    pieces.append(odd_primes)
    return np.concatenate(pieces) if len(pieces) > 1 else pieces[0]


def higher_prime_powers(
    limit: int, base_primes: list[int]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return sorted q=p^a events with a>=2 and their prime metadata."""

    rows: list[tuple[int, int, int, float]] = []
    for prime in base_primes:
        value = prime * prime
        exponent = 2
        log_prime = math.log(prime)
        while value <= limit:
            rows.append((value, prime, exponent, log_prime))
            if value > limit // prime:
                break
            value *= prime
            exponent += 1
    rows.sort(key=lambda row: row[0])
    return (
        np.fromiter((row[0] for row in rows), dtype=np.int64),
        np.fromiter((row[1] for row in rows), dtype=np.int64),
        np.fromiter((row[2] for row in rows), dtype=np.int16),
        np.fromiter((row[3] for row in rows), dtype=np.float64),
    )


def trivial_f_array(y: np.ndarray) -> np.ndarray:
    """Evaluate F(y)=2*atanh(y)-2*y+y*log(1-y^2) stably."""

    output = np.empty_like(y)
    small = y < 0.01
    if np.any(small):
        values = y[small]
        square = values * values
        power = values * square
        total = np.zeros_like(values)
        for layer in range(1, 7):
            total -= power / (layer * (2 * layer + 1))
            power *= square
        output[small] = total
    if np.any(~small):
        values = y[~small]
        output[~small] = (
            (1.0 + values) * np.log1p(values)
            + (values - 1.0) * np.log1p(-values)
            - 2.0 * values
        )
    return output


def row_at(
    index: int,
    q: np.ndarray,
    prime: np.ndarray,
    exponent: np.ndarray,
    t: np.ndarray,
    phi: np.ndarray,
    phi_prime_left: np.ndarray,
    radius: np.ndarray,
    jump: np.ndarray,
    m_left: np.ndarray,
) -> dict[str, Any]:
    """Format the full state at one event."""

    q_value = int(q[index])
    q_quarter = math.sqrt(math.sqrt(q_value))
    phi_prime_right = float(phi_prime_left[index] - jump[index])
    h_right = phi_prime_right * phi_prime_right - float(radius[index]) ** 2
    return {
        "q": q_value,
        "prime": int(prime[index]),
        "exponent": int(exponent[index]),
        "t": float(t[index]),
        "Phi": float(phi[index]),
        "Phi_prime_left": float(phi_prime_left[index]),
        "R": float(radius[index]),
        "J": float(jump[index]),
        "m_left": float(m_left[index]),
        "m_right": float(m_left[index] - jump[index]),
        "ratio_m_left_over_J": float(m_left[index] / jump[index]),
        "m_left_over_q_quarter": float(m_left[index] / q_quarter),
        "m_right_over_q_quarter": float((m_left[index] - jump[index]) / q_quarter),
        "m_left_over_R": float(m_left[index] / radius[index]),
        "abs_Phi_prime_over_R": float(abs(phi_prime_left[index]) / radius[index]),
        "H_right": h_right,
        "minus_H_right": -h_right,
        "minus_H_right_over_sqrt_q": -h_right / math.sqrt(q_value),
    }


def new_bin(key: int, left: int, right: int) -> dict[str, Any]:
    """Create an online range-aggregation record."""

    return {
        "key": key,
        "q_left": left,
        "q_right": right,
        "event_count": 0,
        "prime_event_count": 0,
        "higher_power_event_count": 0,
        "minimum_ratio": None,
        "minimum_ratio_row": None,
        "minimum_phi": None,
        "minimum_phi_row": None,
        "minimum_m_left_over_q_quarter": None,
        "minimum_m_left_over_q_quarter_row": None,
        "minimum_m_right": None,
        "minimum_m_right_row": None,
        "minimum_m_right_over_q_quarter": None,
        "minimum_m_right_over_q_quarter_row": None,
        "minimum_m_left_over_R": None,
        "minimum_m_left_over_R_row": None,
        "maximum_abs_phi_prime_over_R": None,
        "maximum_abs_phi_prime_over_R_row": None,
        "maximum_H_right": None,
        "maximum_H_right_row": None,
        "minimum_H_right": None,
        "minimum_H_right_row": None,
        "minimum_minus_H_right_over_sqrt_q": None,
        "minimum_minus_H_right_over_sqrt_q_row": None,
        "minimum_J": None,
        "minimum_J_row": None,
        "maximum_J": None,
        "maximum_J_row": None,
        "transition_count": 0,
        "sum_smooth_gain": 0.0,
        "sum_withdrawal_on_transitions": 0.0,
        "sum_net_reserve_change": 0.0,
        "minimum_smooth_gain": None,
        "minimum_smooth_gain_row": None,
        "minimum_smooth_gain_over_J": None,
        "minimum_smooth_gain_over_J_row": None,
    }


def update_minimum(
    record: dict[str, Any], name: str, value: float, row: dict[str, Any]
) -> None:
    """Update an online minimum and its witness row."""

    current = record[name]
    if current is None or value < current:
        record[name] = float(value)
        record[f"{name}_row"] = row


def update_maximum(
    record: dict[str, Any], name: str, value: float, row: dict[str, Any]
) -> None:
    """Update an online maximum and its witness row."""

    current = record[name]
    if current is None or value > current:
        record[name] = float(value)
        record[f"{name}_row"] = row


def update_range_bins(
    bins: dict[int, dict[str, Any]],
    bin_keys: np.ndarray,
    q: np.ndarray,
    prime: np.ndarray,
    exponent: np.ndarray,
    t: np.ndarray,
    phi: np.ndarray,
    phi_prime_left: np.ndarray,
    radius: np.ndarray,
    jump: np.ndarray,
    m_left: np.ndarray,
    smooth_gain: np.ndarray,
    base: int,
) -> None:
    """Aggregate event diagnostics into integer logarithmic bins."""

    for key_value in np.unique(bin_keys):
        key = int(key_value)
        positions = np.flatnonzero(bin_keys == key)
        if base == 10:
            left = 10**key
            right = 10 ** (key + 1) - 1
        else:
            left = int(math.floor(math.exp(key)))
            right = int(math.ceil(math.exp(key + 1))) - 1
        record = bins.setdefault(key, new_bin(key, left, right))
        record["event_count"] += int(positions.size)
        local_exponents = exponent[positions]
        record["prime_event_count"] += int(np.count_nonzero(local_exponents == 1))
        record["higher_power_event_count"] += int(
            np.count_nonzero(local_exponents > 1)
        )

        ratio = m_left[positions] / jump[positions]
        phi_local = phi[positions]
        q_quarter = np.sqrt(np.sqrt(q[positions].astype(np.float64)))
        normalized = m_left[positions] / q_quarter
        m_right = m_left[positions] - jump[positions]
        normalized_right = m_right / q_quarter
        m_over_radius = m_left[positions] / radius[positions]
        derivative_ratio = np.abs(phi_prime_left[positions]) / radius[positions]
        phi_prime_right = phi_prime_left[positions] - jump[positions]
        h_right = phi_prime_right * phi_prime_right - radius[positions] ** 2
        normalized_negative_h = -h_right / np.sqrt(q[positions].astype(np.float64))

        local_metrics = [
            ("minimum_ratio", ratio, np.argmin),
            ("minimum_phi", phi_local, np.argmin),
            ("minimum_m_left_over_q_quarter", normalized, np.argmin),
            ("minimum_m_right", m_right, np.argmin),
            ("minimum_m_right_over_q_quarter", normalized_right, np.argmin),
            ("minimum_m_left_over_R", m_over_radius, np.argmin),
            ("maximum_abs_phi_prime_over_R", derivative_ratio, np.argmax),
            ("maximum_H_right", h_right, np.argmax),
            ("minimum_H_right", h_right, np.argmin),
            (
                "minimum_minus_H_right_over_sqrt_q",
                normalized_negative_h,
                np.argmin,
            ),
            ("minimum_J", jump[positions], np.argmin),
            ("maximum_J", jump[positions], np.argmax),
        ]
        for name, values, selector in local_metrics:
            local_position = int(selector(values))
            index = int(positions[local_position])
            row = row_at(
                index, q, prime, exponent, t, phi, phi_prime_left, radius, jump, m_left
            )
            value = float(values[local_position])
            if name.startswith("maximum"):
                update_maximum(record, name, value, row)
            else:
                update_minimum(record, name, value, row)

        transition_positions = positions[np.isfinite(smooth_gain[positions])]
        if transition_positions.size:
            gains = smooth_gain[transition_positions]
            withdrawals = jump[transition_positions]
            record["transition_count"] += int(transition_positions.size)
            record["sum_smooth_gain"] += float(np.sum(gains, dtype=np.float64))
            record["sum_withdrawal_on_transitions"] += float(
                np.sum(withdrawals, dtype=np.float64)
            )
            record["sum_net_reserve_change"] += float(
                np.sum(gains - withdrawals, dtype=np.float64)
            )
            gain_position = int(np.argmin(gains))
            gain_ratio = gains / withdrawals
            ratio_position = int(np.argmin(gain_ratio))
            for name, values, local_position in [
                ("minimum_smooth_gain", gains, gain_position),
                ("minimum_smooth_gain_over_J", gain_ratio, ratio_position),
            ]:
                index = int(transition_positions[local_position])
                row = row_at(
                    index,
                    q,
                    prime,
                    exponent,
                    t,
                    phi,
                    phi_prime_left,
                    radius,
                    jump,
                    m_left,
                )
                row["smooth_gain_from_previous_event_right"] = float(
                    smooth_gain[index]
                )
                row["smooth_gain_over_J"] = float(smooth_gain[index] / jump[index])
                update_minimum(record, name, float(values[local_position]), row)


def finalize_bins(bins: dict[int, dict[str, Any]]) -> list[dict[str, Any]]:
    """Convert online bin records to sorted JSON rows."""

    rows: list[dict[str, Any]] = []
    for key in sorted(bins):
        row = bins[key]
        withdrawals = row["sum_withdrawal_on_transitions"]
        row["smooth_gain_over_withdrawal_sum"] = (
            row["sum_smooth_gain"] / withdrawals if withdrawals else None
        )
        row["balance_identity_residual"] = (
            row["sum_smooth_gain"]
            - withdrawals
            - row["sum_net_reserve_change"]
        )
        rows.append(row)
    return rows


def merge_tightest(
    heap: list[tuple[float, int, dict[str, Any]]],
    rows: list[dict[str, Any]],
    keep: int,
) -> None:
    """Retain the globally smallest reserve-to-kick ratios."""

    for row in rows:
        ratio = float(row["ratio_m_left_over_J"])
        item = (-ratio, int(row["q"]), row)
        if len(heap) < keep:
            heapq.heappush(heap, item)
        elif ratio < -heap[0][0]:
            heapq.heapreplace(heap, item)


def new_monotonicity_record(definition: str) -> dict[str, Any]:
    """Create an online monotonicity diagnostic record."""

    return {
        "definition": definition,
        "previous_value": None,
        "transition_count": 0,
        "increase_count": 0,
        "decrease_count": 0,
        "equal_count": 0,
        "largest_increase": None,
        "largest_increase_row": None,
        "largest_drop": None,
        "largest_drop_row": None,
    }


def update_monotonicity_record(
    record: dict[str, Any],
    values: np.ndarray,
    q: np.ndarray,
    prime: np.ndarray,
    exponent: np.ndarray,
    t: np.ndarray,
    phi: np.ndarray,
    phi_prime_left: np.ndarray,
    radius: np.ndarray,
    jump: np.ndarray,
    m_left: np.ndarray,
) -> None:
    """Update raw binary64 event-to-event monotonicity counts and witnesses."""

    differences = np.empty_like(values)
    previous = record["previous_value"]
    if previous is None:
        differences[0] = np.nan
    else:
        differences[0] = values[0] - float(previous)
    if values.size > 1:
        differences[1:] = values[1:] - values[:-1]
    finite = np.isfinite(differences)
    local = differences[finite]
    record["transition_count"] += int(local.size)
    record["increase_count"] += int(np.count_nonzero(local > 0.0))
    record["decrease_count"] += int(np.count_nonzero(local < 0.0))
    record["equal_count"] += int(np.count_nonzero(local == 0.0))

    if local.size:
        positions = np.flatnonzero(finite)
        increase_position = int(positions[int(np.argmax(local))])
        drop_position = int(positions[int(np.argmin(local))])
        largest_increase = float(differences[increase_position])
        largest_drop = float(differences[drop_position])
        if (
            record["largest_increase"] is None
            or largest_increase > record["largest_increase"]
        ):
            record["largest_increase"] = largest_increase
            record["largest_increase_row"] = row_at(
                increase_position,
                q,
                prime,
                exponent,
                t,
                phi,
                phi_prime_left,
                radius,
                jump,
                m_left,
            )
        if record["largest_drop"] is None or largest_drop < record["largest_drop"]:
            record["largest_drop"] = largest_drop
            record["largest_drop_row"] = row_at(
                drop_position,
                q,
                prime,
                exponent,
                t,
                phi,
                phi_prime_left,
                radius,
                jump,
                m_left,
            )
    record["previous_value"] = float(values[-1])


def finalize_monotonicity_record(record: dict[str, Any]) -> dict[str, Any]:
    """Remove online state and classify a monotonicity diagnostic."""

    output = {key: value for key, value in record.items() if key != "previous_value"}
    output["nondecreasing_on_scanned_events"] = output["decrease_count"] == 0
    output["nonincreasing_on_scanned_events"] = output["increase_count"] == 0
    return output


def scan(args: argparse.Namespace) -> dict[str, Any]:
    """Run the segmented event scan and return a JSON-serializable report."""

    started = time.perf_counter()
    base_primes = simple_primes(math.isqrt(args.limit))
    higher_q, higher_prime, higher_exponent, higher_log = higher_prime_powers(
        args.limit, base_primes
    )

    decade_bins: dict[int, dict[str, Any]] = {}
    natural_log_bins: dict[int, dict[str, Any]] = {}
    tightest_heap: list[tuple[float, int, dict[str, Any]]] = []
    envelope_powers = [0.0, 0.5, 1.0, 2.0]
    envelope_minima: dict[float, tuple[float, dict[str, Any]]] = {}
    round_constants = [0.05, 0.075, 0.1, 0.125, 0.15]
    round_constant_violations = {constant: 0 for constant in round_constants}
    monotonicity = {
        "post_m": new_monotonicity_record("m(log(q)+)"),
        "post_m_over_q_quarter": new_monotonicity_record(
            "m(log(q)+)/q^(1/4)"
        ),
        "pre_m_over_q_quarter": new_monotonicity_record(
            "m(log(q)-)/q^(1/4)"
        ),
        "pre_m_over_J": new_monotonicity_record("m(log(q)-)/J_q"),
        "minus_H_right": new_monotonicity_record("-H(log(q)+)"),
        "minus_H_right_over_sqrt_q": new_monotonicity_record(
            "-H(log(q)+)/sqrt(q)"
        ),
    }

    global_record = new_bin(-1, 2, args.limit)
    prime_count = 0
    event_count = 0
    psi_total = 0.0
    s_lambda_total = 0.0
    previous_m_right: float | None = None
    first_m_right: float | None = None
    last_m_right: float | None = None
    total_smooth_gain = 0.0
    total_transition_withdrawal = 0.0
    max_continuity_error = 0.0
    max_jump_error = 0.0
    max_full_coordinate_b_error = 0.0
    max_full_coordinate_derivative_error = 0.0
    max_bilinear_h_error = 0.0
    previous_h_right: float | None = None
    maximum_open_cell_h_change = -math.inf
    open_cell_h_increase_count_above_1e_minus_8 = 0
    open_cell_transition_count = 0
    min_radius_squared = math.inf
    nonpositive_phi_count = 0
    nonpositive_m_right_count = 0
    processed_high_power_count = 0

    low = 2
    segment_index = 0
    while low <= args.limit:
        high = min(args.limit, low + args.segment_span - 1)
        primes = prime_array_in_segment(low, high, base_primes)
        prime_count += int(primes.size)

        hp_start = int(np.searchsorted(higher_q, low, side="left"))
        hp_stop = int(np.searchsorted(higher_q, high, side="right"))
        local_higher_count = hp_stop - hp_start
        processed_high_power_count += local_higher_count

        prime_logs = np.log(primes.astype(np.float64))
        event_q = np.concatenate((primes, higher_q[hp_start:hp_stop]))
        event_prime = np.concatenate((primes, higher_prime[hp_start:hp_stop]))
        event_exponent = np.concatenate(
            (
                np.ones(primes.size, dtype=np.int16),
                higher_exponent[hp_start:hp_stop],
            )
        )
        event_lambda = np.concatenate((prime_logs, higher_log[hp_start:hp_stop]))
        order = np.argsort(event_q, kind="stable")
        event_q = event_q[order]
        event_prime = event_prime[order]
        event_exponent = event_exponent[order]
        event_lambda = event_lambda[order]
        local_count = int(event_q.size)
        if not local_count:
            low = high + 1
            continue

        q_float = event_q.astype(np.float64)
        reciprocal_weight = event_lambda / q_float
        psi_cumulative = np.cumsum(event_lambda, dtype=np.float64)
        s_cumulative = np.cumsum(reciprocal_weight, dtype=np.float64)
        psi_left = psi_total + psi_cumulative - event_lambda
        s_left = s_lambda_total + s_cumulative - reciprocal_weight
        psi_right = psi_total + psi_cumulative
        s_right = s_lambda_total + s_cumulative

        t = np.log(q_float)
        root_q = np.sqrt(q_float)
        y = 1.0 / q_float
        f_value = trivial_f_array(y)
        f_prime = np.log1p(-(y * y))
        trivial_b = 0.5 * root_q * f_value
        trivial_derivative = 0.25 * root_q * f_value - 0.5 * root_q * y * f_prime

        u_left = s_left - t + EULER_GAMMA
        v_left = psi_left - q_float + LOG_TWO_PI
        b_left = root_q * u_left - v_left / root_q - trivial_b
        b_prime_left = 0.5 * (root_q * u_left + v_left / root_q) - trivial_derivative

        u_right = s_right - t + EULER_GAMMA
        v_right = psi_right - q_float + LOG_TWO_PI
        b_right = root_q * u_right - v_right / root_q - trivial_b
        b_prime_right = (
            0.5 * (root_q * u_right + v_right / root_q) - trivial_derivative
        )
        max_continuity_error = max(
            max_continuity_error, float(np.max(np.abs(b_right - b_left)))
        )

        jump = event_lambda / root_q
        max_jump_error = max(
            max_jump_error,
            float(np.max(np.abs((b_prime_right - b_prime_left) - jump))),
        )
        b_event = 0.5 * (b_left + b_right)
        phi = WEIL_C - b_event
        phi_prime_left = -b_prime_left
        g = root_q * (1.0 - y**3 / (1.0 - y * y))
        k_value = g - 0.25 * WEIL_C
        radius_squared = 0.25 * phi * phi + 2.0 * k_value * phi
        min_radius_squared = min(min_radius_squared, float(np.min(radius_squared)))
        nonpositive_phi_count += int(np.count_nonzero(phi <= 0.0))
        radius = np.sqrt(np.maximum(radius_squared, 0.0))
        m_left = radius + phi_prime_left
        m_right = m_left - jump
        nonpositive_m_right_count += int(np.count_nonzero(m_right <= 0.0))

        # Full-B characteristic coordinates, including the trivial term.
        # The u correction is evaluated through F to avoid cancellation in
        # y-atanh(y) for large q.
        u_full_right = u_right + 0.5 * y * f_prime - 0.5 * f_value
        v_full_right = v_right + 0.5 * f_prime
        coordinate_b = root_q * u_full_right - v_full_right / root_q
        coordinate_b_prime = 0.5 * (
            root_q * u_full_right + v_full_right / root_q
        )
        max_full_coordinate_b_error = max(
            max_full_coordinate_b_error,
            float(np.max(np.abs(coordinate_b - b_right))),
        )
        max_full_coordinate_derivative_error = max(
            max_full_coordinate_derivative_error,
            float(np.max(np.abs(coordinate_b_prime - b_prime_right))),
        )

        phi_prime_right = phi_prime_left - jump
        h_right = phi_prime_right * phi_prime_right - radius_squared
        h_left = phi_prime_left * phi_prime_left - radius_squared
        reciprocal_x2_minus_one = 1.0 / (q_float * q_float - 1.0)
        tiny_g_correction = reciprocal_x2_minus_one / q_float
        h_small = reciprocal_x2_minus_one / root_q
        h_bilinear = (
            u_full_right * v_full_right
            + 2.0 * (q_float - reciprocal_x2_minus_one) * u_full_right
            - 2.0 * (1.0 - tiny_g_correction) * v_full_right
            + 0.25 * WEIL_C * WEIL_C
            - 2.0 * WEIL_C * (root_q - h_small)
        )
        max_bilinear_h_error = max(
            max_bilinear_h_error,
            float(np.max(np.abs(h_bilinear - h_right))),
        )

        open_cell_change = np.empty(local_count, dtype=np.float64)
        if previous_h_right is None:
            open_cell_change[0] = np.nan
        else:
            open_cell_change[0] = h_left[0] - previous_h_right
        if local_count > 1:
            open_cell_change[1:] = h_left[1:] - h_right[:-1]
        finite_open_change = np.isfinite(open_cell_change)
        if np.any(finite_open_change):
            local_open_change = open_cell_change[finite_open_change]
            maximum_open_cell_h_change = max(
                maximum_open_cell_h_change, float(np.max(local_open_change))
            )
            open_cell_h_increase_count_above_1e_minus_8 += int(
                np.count_nonzero(local_open_change > 1e-8)
            )
            open_cell_transition_count += int(local_open_change.size)
        previous_h_right = float(h_right[-1])

        smooth_gain = np.full(local_count, np.nan, dtype=np.float64)
        if previous_m_right is not None:
            smooth_gain[0] = m_left[0] - previous_m_right
        if local_count > 1:
            smooth_gain[1:] = m_left[1:] - m_right[:-1]
        if first_m_right is None:
            first_m_right = float(m_right[0])
        previous_m_right = float(m_right[-1])
        last_m_right = previous_m_right
        finite_gain = np.isfinite(smooth_gain)
        total_smooth_gain += float(np.sum(smooth_gain[finite_gain], dtype=np.float64))
        total_transition_withdrawal += float(
            np.sum(jump[finite_gain], dtype=np.float64)
        )

        decade_keys = np.floor(np.log10(q_float)).astype(np.int16)
        natural_keys = np.floor(t).astype(np.int16)
        update_range_bins(
            decade_bins,
            decade_keys,
            event_q,
            event_prime,
            event_exponent,
            t,
            phi,
            phi_prime_left,
            radius,
            jump,
            m_left,
            smooth_gain,
            10,
        )
        update_range_bins(
            natural_log_bins,
            natural_keys,
            event_q,
            event_prime,
            event_exponent,
            t,
            phi,
            phi_prime_left,
            radius,
            jump,
            m_left,
            smooth_gain,
            math.e,
        )

        all_positions = np.arange(local_count)
        ratio = m_left / jump
        q_quarter = np.sqrt(root_q)
        metrics = [
            ("minimum_ratio", ratio, np.argmin),
            ("minimum_phi", phi, np.argmin),
            ("minimum_m_left_over_q_quarter", m_left / q_quarter, np.argmin),
            ("minimum_m_right", m_right, np.argmin),
            (
                "minimum_m_right_over_q_quarter",
                m_right / q_quarter,
                np.argmin,
            ),
            ("minimum_m_left_over_R", m_left / radius, np.argmin),
            ("maximum_abs_phi_prime_over_R", np.abs(phi_prime_left) / radius, np.argmax),
            ("maximum_H_right", h_right, np.argmax),
            ("minimum_H_right", h_right, np.argmin),
            (
                "minimum_minus_H_right_over_sqrt_q",
                -h_right / root_q,
                np.argmin,
            ),
            ("minimum_J", jump, np.argmin),
            ("maximum_J", jump, np.argmax),
        ]
        for name, values, selector in metrics:
            index = int(selector(values))
            row = row_at(
                index,
                event_q,
                event_prime,
                event_exponent,
                t,
                phi,
                phi_prime_left,
                radius,
                jump,
                m_left,
            )
            value = float(values[index])
            if name.startswith("maximum"):
                update_maximum(global_record, name, value, row)
            else:
                update_minimum(global_record, name, value, row)

        if np.any(finite_gain):
            finite_positions = all_positions[finite_gain]
            local_gains = smooth_gain[finite_gain]
            gain_index = int(finite_positions[int(np.argmin(local_gains))])
            gain_ratios = local_gains / jump[finite_gain]
            gain_ratio_index = int(finite_positions[int(np.argmin(gain_ratios))])
            for name, values, index in [
                ("minimum_smooth_gain", smooth_gain, gain_index),
                ("minimum_smooth_gain_over_J", smooth_gain / jump, gain_ratio_index),
            ]:
                row = row_at(
                    index,
                    event_q,
                    event_prime,
                    event_exponent,
                    t,
                    phi,
                    phi_prime_left,
                    radius,
                    jump,
                    m_left,
                )
                row["smooth_gain_from_previous_event_right"] = float(smooth_gain[index])
                row["smooth_gain_over_J"] = float(smooth_gain[index] / jump[index])
                update_minimum(global_record, name, float(values[index]), row)

        monotonic_values = {
            "post_m": m_right,
            "post_m_over_q_quarter": m_right / q_quarter,
            "pre_m_over_q_quarter": m_left / q_quarter,
            "pre_m_over_J": ratio,
            "minus_H_right": -h_right,
            "minus_H_right_over_sqrt_q": -h_right / root_q,
        }
        for name, values in monotonic_values.items():
            update_monotonicity_record(
                monotonicity[name],
                values,
                event_q,
                event_prime,
                event_exponent,
                t,
                phi,
                phi_prime_left,
                radius,
                jump,
                m_left,
            )

        local_keep = min(args.keep_tightest, local_count)
        tight_positions = np.argpartition(ratio, local_keep - 1)[:local_keep]
        merge_tightest(
            tightest_heap,
            [
                row_at(
                    int(index),
                    event_q,
                    event_prime,
                    event_exponent,
                    t,
                    phi,
                    phi_prime_left,
                    radius,
                    jump,
                    m_left,
                )
                for index in tight_positions
            ],
            args.keep_tightest,
        )

        for power in envelope_powers:
            normalized = m_left * np.power(1.0 + t, power) / q_quarter
            index = int(np.argmin(normalized))
            value = float(normalized[index])
            previous = envelope_minima.get(power)
            if previous is None or value < previous[0]:
                envelope_minima[power] = (
                    value,
                    row_at(
                        index,
                        event_q,
                        event_prime,
                        event_exponent,
                        t,
                        phi,
                        phi_prime_left,
                        radius,
                        jump,
                        m_left,
                    ),
                )
        for constant in round_constants:
            round_constant_violations[constant] += int(
                np.count_nonzero(m_left < constant * q_quarter)
            )

        event_count += local_count
        global_record["event_count"] += local_count
        global_record["prime_event_count"] += int(primes.size)
        global_record["higher_power_event_count"] += local_higher_count
        psi_total = float(psi_right[-1])
        s_lambda_total = float(s_right[-1])

        segment_index += 1
        if args.progress_every and (
            segment_index % args.progress_every == 0 or high == args.limit
        ):
            elapsed = time.perf_counter() - started
            print(
                f"segment {segment_index}: q<={high:,}, events={event_count:,}, "
                f"min_ratio={global_record['minimum_ratio']:.12g}, "
                f"elapsed={elapsed:.1f}s",
                flush=True,
            )
        low = high + 1

    expected_balance = (
        (last_m_right - first_m_right)
        if first_m_right is not None and last_m_right is not None
        else math.nan
    )
    observed_balance = total_smooth_gain - total_transition_withdrawal
    decade_rows = finalize_bins(decade_bins)
    natural_rows = finalize_bins(natural_log_bins)

    tail_rows: list[dict[str, Any]] = []
    running_ratio: tuple[float, dict[str, Any]] | None = None
    running_phi: tuple[float, dict[str, Any]] | None = None
    running_normalized: tuple[float, dict[str, Any]] | None = None
    for row in reversed(decade_rows):
        candidates = [
            ("minimum_ratio", "minimum_ratio_row", running_ratio),
            ("minimum_phi", "minimum_phi_row", running_phi),
            (
                "minimum_m_left_over_q_quarter",
                "minimum_m_left_over_q_quarter_row",
                running_normalized,
            ),
        ]
        updated: list[tuple[float, dict[str, Any]]] = []
        for value_key, row_key, running in candidates:
            candidate = (float(row[value_key]), row[row_key])
            if running is None or candidate[0] < running[0]:
                running = candidate
            updated.append(running)
        running_ratio, running_phi, running_normalized = updated
        tail_rows.append(
            {
                "q_at_least": int(row["q_left"]),
                "minimum_ratio_m_left_over_J": running_ratio[0],
                "minimum_ratio_row": running_ratio[1],
                "minimum_Phi": running_phi[0],
                "minimum_Phi_row": running_phi[1],
                "minimum_m_left_over_q_quarter": running_normalized[0],
                "minimum_m_left_over_q_quarter_row": running_normalized[1],
            }
        )
    tail_rows.reverse()

    envelope_rows = []
    for power in envelope_powers:
        value, witness = envelope_minima[power]
        envelope_rows.append(
            {
                "beta": power,
                "tested_form": "m_left >= c*q^(1/4)/(1+log(q))^beta",
                "largest_observed_finite_scan_constant_c": value,
                "witness": witness,
            }
        )

    tightest_rows = sorted(
        (item[2] for item in tightest_heap),
        key=lambda row: (row["ratio_m_left_over_J"], row["q"]),
    )
    checks = {
        "all_prime_power_events_merged": processed_high_power_count == higher_q.size,
        "event_count_decomposes": event_count == prime_count + higher_q.size,
        "Phi_positive_at_all_events": nonpositive_phi_count == 0,
        "post_kick_reserve_positive_at_all_events": nonpositive_m_right_count == 0,
        "radius_squared_positive_at_all_events": min_radius_squared > 0.0,
        "smooth_gain_positive_after_first_event": (
            global_record["minimum_smooth_gain"] is not None
            and global_record["minimum_smooth_gain"] > -1e-9
        ),
        "continuity_residual_below_1e_minus_7": max_continuity_error < 1e-7,
        "jump_residual_below_1e_minus_7": max_jump_error < 1e-7,
        "full_coordinate_reconstruction_below_1e_minus_7": max(
            max_full_coordinate_b_error, max_full_coordinate_derivative_error
        )
        < 1e-7,
        "bilinear_H_factorization_residual_below_1e_minus_3": (
            max_bilinear_h_error < 1e-3
        ),
        "raw_open_cell_H_subtraction_residual_below_1e_minus_3": (
            maximum_open_cell_h_change < 1e-3
        ),
        "reserve_balance_residual_below_1e_minus_7": abs(
            observed_balance - expected_balance
        )
        < 1e-7,
    }
    checks["all_pass"] = all(checks.values())

    return {
        "status": (
            "Exploratory binary64 segmented scan. This is finite numerical "
            "evidence, not a directed-rounding certificate, a universal lower "
            "envelope, or a proof of RH."
        ),
        "configuration": {
            "limit": args.limit,
            "segment_span": args.segment_span,
            "keep_tightest": args.keep_tightest,
            "base_prime_limit": math.isqrt(args.limit),
        },
        "exact_quantities_evaluated": {
            "Phi": "Phi=C-B, C=2+gamma_E-log(4*pi)",
            "R": "R=sqrt(Phi^2/4+2*k*Phi), k=g-C/4",
            "g": "g=sqrt(q)*(1-q^(-3)/(1-q^(-2)))",
            "m": "m=R+Phi'",
            "event_target": "m(log(q)-)>J_q, J_q=Lambda(q)/sqrt(q)",
            "event_update": "m(log(q)+)=m(log(q)-)-J_q",
            "smooth_balance": (
                "gain_i=m(log(q_i)-)-m(log(q_(i-1))+), positive while the "
                "first-crossing cone hypotheses hold"
            ),
        },
        "counts": {
            "prime_count": prime_count,
            "higher_prime_power_count": int(higher_q.size),
            "all_prime_power_event_count": event_count,
            "nonpositive_Phi_count": nonpositive_phi_count,
            "nonpositive_post_kick_reserve_count": nonpositive_m_right_count,
        },
        "global_extrema": global_record,
        "tightest_event_ratios": tightest_rows,
        "decade_diagnostics": decade_rows,
        "natural_log_bin_diagnostics": natural_rows,
        "tail_minima_at_decade_thresholds": tail_rows,
        "candidate_q_quarter_envelopes": envelope_rows,
        "round_constant_tests": [
            {
                "candidate": f"m_left >= {constant:g}*q^(1/4)",
                "violation_count": round_constant_violations[constant],
                "passes_finite_scan": round_constant_violations[constant] == 0,
            }
            for constant in round_constants
        ],
        "event_to_event_balance": {
            "first_event_post_kick_m": first_m_right,
            "last_event_post_kick_m": last_m_right,
            "sum_smooth_gains_after_first_event": total_smooth_gain,
            "sum_withdrawals_after_first_event": total_transition_withdrawal,
            "gain_minus_withdrawal": observed_balance,
            "last_minus_first_post_kick_m": expected_balance,
            "telescoping_residual": observed_balance - expected_balance,
            "interpretation": (
                "Both cumulative terms are much larger than their difference; "
                "a proof needs correlated cancellation, not separate crude bounds."
            ),
        },
        "bilinear_post_event_energy": {
            "coordinates": (
                "u=S_Lambda-log(x)+gamma_E+1/x-atanh(1/x), "
                "v=psi-x+log(2*pi)+log(1-x^(-2))/2"
            ),
            "factorization": (
                "H_+=u*v+2*(x-1/(x^2-1))*u"
                "-2*(1-1/(x*(x^2-1)))*v+C^2/4"
                "-2*C*(sqrt(x)-h(x))"
            ),
            "h": "h(x)=1/(sqrt(x)*(x^2-1))",
            "derivation": (
                "B=sqrt(x)*u-v/sqrt(x), B'=(sqrt(x)*u+v/sqrt(x))/2; "
                "substitution into H=Phi'^2-Phi^2/4-2*k*Phi gives the formula."
            ),
            "maximum_absolute_binary64_factorization_residual": max_bilinear_h_error,
        },
        "event_sequence_monotonicity": {
            name: finalize_monotonicity_record(record)
            for name, record in monotonicity.items()
        },
        "local_open_cell_barrier": {
            "exact_law": "H'=-2*k'(t)*Phi(t)<0 while Phi>0",
            "tested_transition_count": open_cell_transition_count,
            "maximum_observed_H_left_next_minus_H_right_previous": (
                maximum_open_cell_h_change
            ),
            "increase_count_above_1e_minus_8": (
                open_cell_h_increase_count_above_1e_minus_8
            ),
            "interpretation": (
                "The exact derivative law proves the local barrier under Phi>0. "
                "Direct subtraction of nearly equal binary64 endpoint energies "
                "is ill-conditioned on very short cells, so small positive raw "
                "residuals are reported as roundoff diagnostics rather than sign "
                "violations. The unresolved step remains preservation across the "
                "event kick."
            ),
        },
        "numerical_validation": {
            "max_abs_B_right_minus_B_left": max_continuity_error,
            "max_abs_derivative_jump_minus_J": max_jump_error,
            "max_abs_full_coordinate_B_minus_direct_B": max_full_coordinate_b_error,
            "max_abs_full_coordinate_B_prime_minus_direct_B_prime": (
                max_full_coordinate_derivative_error
            ),
            "max_abs_bilinear_H_minus_energy_H": max_bilinear_h_error,
            "minimum_radius_squared": min_radius_squared,
        },
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
    parser.add_argument("--limit", type=int, default=100_000_000)
    parser.add_argument("--segment-span", type=int, default=5_000_000)
    parser.add_argument("--keep-tightest", type=int, default=25)
    parser.add_argument("--progress-every", type=int, default=1)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()
    if args.limit < 100:
        raise ValueError("limit must be at least 100")
    if args.segment_span < 10_000:
        raise ValueError("segment-span must be at least 10,000")
    if args.keep_tightest < 1:
        raise ValueError("keep-tightest must be positive")

    report = scan(args)
    if not report["checks"]["all_pass"]:
        raise AssertionError(f"exploratory consistency checks failed: {report['checks']}")
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.report}")
    print(f"all checks pass: {report['checks']['all_pass']}")
    print(f"events: {report['counts']['all_prime_power_event_count']:,}")
    print(
        "minimum m_left/J: "
        f"{report['global_extrema']['minimum_ratio']:.15g} at "
        f"q={report['global_extrema']['minimum_ratio_row']['q']:,}"
    )
    print(
        "minimum m_left/q^(1/4): "
        f"{report['global_extrema']['minimum_m_left_over_q_quarter']:.15g} at "
        f"q={report['global_extrema']['minimum_m_left_over_q_quarter_row']['q']:,}"
    )


if __name__ == "__main__":
    main()
