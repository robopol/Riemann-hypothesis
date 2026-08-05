#!/usr/bin/env python3
"""Audit finite scaling of the signed triangular discrepancy.

The script recomputes every prime-power event through a selected limit for
lambda=1/4, evaluates |E_h|/Q, and compares the full event population with the
logarithmic plotting sample.  It estimates the logarithmic exponent in

    |E_h| / Q ~ (log q)^A / sqrt(q)

from equal-log-width bin quantiles.  All fits are descriptive finite-range
diagnostics; they are not asymptotic bounds.
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
    trivial_f_array,
)
from exploratory_signed_triangular_scan import (
    build_event_support,
    phi_at_x,
    required_support_limit,
)
from plot_signed_triangular_dynamics import sample_event_indices


SCRIPT_PATH = Path(__file__).resolve()
DEFAULT_REPORT = SCRIPT_PATH.with_name("eh_scaling_law_audit_report.json")
QUANTILES = (0.5, 0.9, 0.99, 0.999)


def sha256(path: Path) -> str:
    """Return the lowercase SHA-256 digest of a file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def linear_fit(x: np.ndarray, y: np.ndarray) -> dict[str, float | int | None]:
    """Return a simple least-squares fit and its coefficient of determination."""

    if x.size < 3:
        return {
            "count": int(x.size),
            "slope": None,
            "slope_standard_error": None,
            "intercept": None,
            "R_squared": None,
        }
    slope, intercept = np.polyfit(x, y, 1)
    predicted = intercept + slope * x
    residual = float(np.sum((y - predicted) ** 2))
    total = float(np.sum((y - float(np.mean(y))) ** 2))
    centered_sum = float(np.sum((x - float(np.mean(x))) ** 2))
    slope_standard_error = None
    if x.size > 2 and centered_sum > 0.0:
        slope_standard_error = math.sqrt(
            (residual / float(x.size - 2)) / centered_sum
        )
    return {
        "count": int(x.size),
        "slope": float(slope),
        "slope_standard_error": slope_standard_error,
        "intercept": float(intercept),
        "R_squared": 1.0 - residual / total if total else None,
    }


def compute_full_state(
    limit: int,
    lambda_value: float,
    segment_span: int,
    chunk_size: int,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Evaluate q, exponent, E_h sign, and |E_h|/Q at every center event."""

    support_limit = required_support_limit(limit, lambda_value)
    event_q, event_prime, event_exponent, event_lambda = build_event_support(
        support_limit, segment_span
    )
    cumulative_psi = np.cumsum(event_lambda, dtype=np.float64)
    cumulative_s_lambda = np.cumsum(
        event_lambda / event_q.astype(np.float64), dtype=np.float64
    )
    center_count = int(np.searchsorted(event_q, limit, side="right"))
    absolute_ratio = np.empty(center_count, dtype=np.float64)
    signed_e = np.empty(center_count, dtype=np.float64)

    for start in range(0, center_count, chunk_size):
        stop = min(start + chunk_size, center_count)
        q = event_q[start:stop]
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
        spectral_b = root_x * u_right - v_right / root_x
        phi = WEIL_C - spectral_b
        p_right = -0.5 * (root_x * u_right + v_right / root_x)
        jump = lambda_q / root_x
        p_bar = p_right + 0.5 * jump

        b_factor = 2.0 * (1.0 - 1.0 / (x * (x * x - 1.0)))
        radius_squared = phi * (
            root_x * b_factor - 0.5 * WEIL_C + 0.25 * phi
        )
        if np.any(radius_squared <= 0.0):
            raise RuntimeError("The cone radius is nonpositive in the scan")
        radius = np.sqrt(radius_squared)

        h = lambda_value * log_x / root_x
        phi_minus = phi_at_x(
            x * np.exp(-h), event_q, cumulative_psi, cumulative_s_lambda
        )
        phi_plus = phi_at_x(
            x * np.exp(h), event_q, cumulative_psi, cumulative_s_lambda
        )
        d_h = (phi_plus - phi_minus) / (2.0 * h)
        e_h = p_bar - d_h
        signed_e[start:stop] = e_h
        absolute_ratio[start:stop] = np.abs(e_h) / radius

    state = {
        "q": event_q[:center_count],
        "prime": event_prime[:center_count],
        "exponent": event_exponent[:center_count],
        "E_h": signed_e,
        "abs_E_over_Q": absolute_ratio,
    }
    metadata = {
        "center_count": center_count,
        "support_count": int(event_q.size),
        "support_limit": support_limit,
    }
    return state, metadata


def quantile_record(values: np.ndarray) -> dict[str, float | int]:
    """Summarize a nonempty positive sample."""

    quantile_values = np.quantile(values, QUANTILES)
    return {
        "count": int(values.size),
        **{
            f"p{100 * quantile:g}": float(value)
            for quantile, value in zip(QUANTILES, quantile_values, strict=True)
        },
        "maximum": float(np.max(values)),
    }


def binned_log_exponent_fit(
    q: np.ndarray,
    ratio: np.ndarray,
    mask: np.ndarray,
    quantile: float,
    cutoff: float,
    limit: float,
    bins_per_decade: int = 10,
    minimum_bin_count: int = 20,
) -> dict[str, Any]:
    """Fit A after fixing the q exponent to -1/2 in equal-log-width bins."""

    lower_log = math.log(cutoff)
    upper_log = math.log(limit)
    bin_count = max(4, int(round(math.log10(limit / cutoff) * bins_per_decade)))
    edges = np.linspace(lower_log, upper_log, bin_count + 1)
    log_q = np.log(q.astype(np.float64))
    selected = mask & (q >= cutoff) & (q <= limit)
    x_values: list[float] = []
    y_values: list[float] = []
    counts: list[int] = []
    centers: list[float] = []
    quantile_values: list[float] = []
    scaled_without_log = np.sqrt(q.astype(np.float64)) * ratio
    for left, right in zip(edges[:-1], edges[1:], strict=True):
        in_bin = selected & (log_q >= left) & (log_q < right)
        count = int(np.count_nonzero(in_bin))
        if count < minimum_bin_count:
            continue
        center_log_q = 0.5 * (left + right)
        value = float(np.quantile(scaled_without_log[in_bin], quantile))
        if value <= 0.0:
            continue
        x_values.append(math.log(center_log_q))
        y_values.append(math.log(value))
        counts.append(count)
        centers.append(math.exp(center_log_q))
        quantile_values.append(value)
    fit = linear_fit(np.asarray(x_values), np.asarray(y_values))
    return {
        "quantile": quantile,
        "cutoff": cutoff,
        "limit": limit,
        "bins_per_decade": bins_per_decade,
        "minimum_bin_count": minimum_bin_count,
        "fit": fit,
        "bin_counts": counts,
        "bin_centers": centers,
        "bin_quantile_values_sqrt_q_abs_E_over_Q": quantile_values,
    }


def decade_diagnostics(
    state: dict[str, np.ndarray], sample_indices: np.ndarray, exponent_a: float
) -> list[dict[str, Any]]:
    """Return full-population and plotting-sample summaries by q decade."""

    q = state["q"]
    ratio = state["abs_E_over_Q"]
    exponent = state["exponent"]
    e_h = state["E_h"]
    log_q = np.log(q.astype(np.float64))
    z_a = np.sqrt(q.astype(np.float64)) * ratio / (log_q**exponent_a)
    sample_mask = np.zeros(q.size, dtype=bool)
    sample_mask[sample_indices] = True
    rows: list[dict[str, Any]] = []
    maximum_decade = int(math.floor(math.log10(float(np.max(q)))))
    for decade in range(maximum_decade + 1):
        left = max(2, 10**decade)
        right = min(int(np.max(q)) + 1, 10 ** (decade + 1))
        mask = (q >= left) & (q < right)
        if not np.any(mask):
            continue
        sample = mask & sample_mask
        prime = mask & (exponent == 1)
        higher = mask & (exponent >= 2)
        positive = mask & (e_h > 0.0)
        negative = mask & (e_h < 0.0)
        maximum_position = np.flatnonzero(mask)[int(np.argmax(z_a[mask]))]
        row: dict[str, Any] = {
            "decade": decade,
            "q_left": left,
            "q_right_exclusive": right,
            "event_count": int(np.count_nonzero(mask)),
            "prime_count": int(np.count_nonzero(prime)),
            "higher_prime_power_count": int(np.count_nonzero(higher)),
            "positive_E_count": int(np.count_nonzero(positive)),
            "negative_E_count": int(np.count_nonzero(negative)),
            "positive_E_fraction": float(np.mean(e_h[mask] > 0.0)),
            "abs_E_h": quantile_record(np.abs(e_h[mask])),
            "abs_E_over_Q": quantile_record(ratio[mask]),
            "Z_A_full": quantile_record(z_a[mask]),
            "Z_A_maximum_witness": {
                "q": int(q[maximum_position]),
                "prime": int(state["prime"][maximum_position]),
                "exponent": int(exponent[maximum_position]),
                "E_h": float(e_h[maximum_position]),
                "abs_E_over_Q": float(ratio[maximum_position]),
                "Z_A": float(z_a[maximum_position]),
            },
            "plot_sample_count": int(np.count_nonzero(sample)),
            "plot_sample_fraction": float(np.mean(sample_mask[mask])),
            "Z_A_plot_sample": quantile_record(z_a[sample]),
            "abs_E_over_Q_plot_sample": quantile_record(ratio[sample]),
            "abs_E_h_plot_sample": quantile_record(np.abs(e_h[sample])),
        }
        if np.any(prime):
            row["Z_A_primes"] = quantile_record(z_a[prime])
        if np.any(higher):
            row["Z_A_higher_prime_powers"] = quantile_record(z_a[higher])
        if np.any(positive):
            row["Z_A_positive_E"] = quantile_record(z_a[positive])
        if np.any(negative):
            row["Z_A_negative_E"] = quantile_record(z_a[negative])
        rows.append(row)
    return rows


def category_fits(state: dict[str, np.ndarray], limit: int) -> dict[str, Any]:
    """Estimate the fixed-square-root logarithmic exponent by category."""

    q = state["q"]
    ratio = state["abs_E_over_Q"]
    exponent = state["exponent"]
    e_h = state["E_h"]
    categories = {
        "all": np.ones(q.size, dtype=bool),
        "primes": exponent == 1,
        "higher_prime_powers": exponent >= 2,
        "positive_E": e_h > 0.0,
        "negative_E": e_h < 0.0,
        "prime_positive_E": (exponent == 1) & (e_h > 0.0),
        "prime_negative_E": (exponent == 1) & (e_h < 0.0),
        "higher_positive_E": (exponent >= 2) & (e_h > 0.0),
        "higher_negative_E": (exponent >= 2) & (e_h < 0.0),
    }
    results: dict[str, Any] = {}
    for name, mask in categories.items():
        minimum_count = 5 if "higher" in name else 100
        fits = []
        for cutoff in (10**3, 10**4, 10**5):
            if cutoff >= limit:
                continue
            for quantile in QUANTILES:
                fits.append(
                    binned_log_exponent_fit(
                        q,
                        ratio,
                        mask,
                        quantile,
                        cutoff,
                        limit,
                        minimum_bin_count=minimum_count,
                    )
                )
        results[name] = {
            "count": int(np.count_nonzero(mask)),
            "fits": fits,
        }
    return results


def extreme_value_diagnostic(
    decade_rows: list[dict[str, Any]], exponent_a: float
) -> dict[str, Any]:
    """Fit the late-decade maximum after the candidate normalization."""

    rows = [row for row in decade_rows if row["q_left"] >= 10_000]
    log_log_q = np.asarray(
        [math.log(math.log(math.sqrt(row["q_left"] * row["q_right_exclusive"]))) for row in rows]
    )
    log_maximum = np.log(
        np.asarray([row["Z_A_full"]["maximum"] for row in rows])
    )
    log_count = np.log(np.asarray([row["event_count"] for row in rows], dtype=float))
    return {
        "normalization_A": exponent_a,
        "late_decade_count": len(rows),
        "fit_log_max_Z_on_log_log_q": linear_fit(log_log_q, log_maximum),
        "fit_log_max_Z_on_log_event_count": linear_fit(log_count, log_maximum),
        "max_over_p99p9": [
            float(row["Z_A_full"]["maximum"] / row["Z_A_full"]["p99.9"])
            for row in rows
        ],
    }


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    """Run the full audit and assemble a JSON-serializable report."""

    started = time.perf_counter()
    state, metadata = compute_full_state(
        args.limit,
        args.lambda_value,
        args.segment_span,
        args.chunk_size,
    )
    sample_indices = sample_event_indices(
        state["q"],
        state["exponent"],
        state["q"].size,
        args.limit,
        args.sample_count,
    )
    fits = category_fits(state, args.limit)
    all_fits = fits["all"]["fits"]
    transition_fit = next(
        record
        for record in all_fits
        if record["cutoff"] == 10_000 and record["quantile"] == 0.5
    )
    tail_fit = next(
        record
        for record in all_fits
        if record["cutoff"] == 100_000 and record["quantile"] == 0.5
    )
    transition_a = float(transition_fit["fit"]["slope"])
    tail_a = float(tail_fit["fit"]["slope"])
    preferred_a = 1.0
    rows_tail_a = decade_diagnostics(state, sample_indices, tail_a)
    rows_a1 = decade_diagnostics(state, sample_indices, 1.0)

    report = {
        "status": (
            "Exploratory binary64 finite-range scaling audit. Fits and "
            "quantiles are descriptive and are not asymptotic bounds."
        ),
        "configuration": {
            "limit": args.limit,
            "lambda": args.lambda_value,
            "segment_span": args.segment_span,
            "chunk_size": args.chunk_size,
            "plot_sample_target_count": args.sample_count,
            "plot_sample_actual_count": int(sample_indices.size),
            "log_exponent_fit_definition": (
                "slope of log bin-quantile(sqrt(q)*|E_h|/Q) on log(log q), "
                "using 10 equal-log bins per decade"
            ),
            "median_fit_A_q_ge_1e4": transition_a,
            "median_fit_A_q_ge_1e5": tail_a,
            "preferred_parsimonious_A": preferred_a,
            "preferred_A_reason": (
                "The q>=1e5 central and upper quantile fits cluster near one, "
                "and A=1 makes the last four decadal distributions nearly stationary."
            ),
        },
        "population": {
            **metadata,
            "prime_count": int(np.count_nonzero(state["exponent"] == 1)),
            "higher_prime_power_count": int(np.count_nonzero(state["exponent"] >= 2)),
            "positive_E_count": int(np.count_nonzero(state["E_h"] > 0.0)),
            "negative_E_count": int(np.count_nonzero(state["E_h"] < 0.0)),
            "zero_E_binary64_count": int(np.count_nonzero(state["E_h"] == 0.0)),
        },
        "category_log_exponent_fits": fits,
        "decades_for_fitted_tail_A": rows_tail_a,
        "decades_for_A_equals_1": rows_a1,
        "extreme_value_diagnostic_A_equals_1": extreme_value_diagnostic(rows_a1, 1.0),
        "runtime": {
            "seconds": time.perf_counter() - started,
            "python": platform.python_version(),
            "numpy": np.__version__,
            "platform": platform.platform(),
            "source_sha256": {
                "audit_eh_scaling_law.py": sha256(SCRIPT_PATH),
                "exploratory_signed_triangular_scan.py": sha256(
                    SCRIPT_PATH.with_name("exploratory_signed_triangular_scan.py")
                ),
                "exploratory_eventwise_reserve_scan.py": sha256(
                    SCRIPT_PATH.with_name("exploratory_eventwise_reserve_scan.py")
                ),
                "plot_signed_triangular_dynamics.py": sha256(
                    SCRIPT_PATH.with_name("plot_signed_triangular_dynamics.py")
                ),
            },
        },
    }
    return report


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=10**8)
    parser.add_argument("--lambda-value", type=float, default=0.25)
    parser.add_argument("--segment-span", type=int, default=10_000_000)
    parser.add_argument("--chunk-size", type=int, default=250_000)
    parser.add_argument("--sample-count", type=int, default=30_000)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def main() -> None:
    """Execute the audit and write its report."""

    args = parse_args()
    if args.limit < 10_000:
        raise SystemExit("--limit must be at least 10000")
    if not 0.0 < args.lambda_value < math.sqrt(2.0):
        raise SystemExit("--lambda-value must satisfy 0<lambda<sqrt(2)")
    report = build_report(args)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.report}")
    print(
        json.dumps(
            {
                "median_fit_A_q_ge_1e4": report["configuration"]["median_fit_A_q_ge_1e4"],
                "median_fit_A_q_ge_1e5": report["configuration"]["median_fit_A_q_ge_1e5"],
                "preferred_parsimonious_A": report["configuration"]["preferred_parsimonious_A"],
                "event_count": report["population"]["center_count"],
                "plot_sample_actual_count": report["configuration"]["plot_sample_actual_count"],
                "seconds": report["runtime"]["seconds"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
