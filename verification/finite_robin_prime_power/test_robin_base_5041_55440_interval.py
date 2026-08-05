from __future__ import annotations

import argparse
import decimal
import hashlib
import json
import math
import platform
import sys
import time
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path

from test_ca_all_profile_interval_certificate import (
    DecimalIntervals,
    harmonic_gamma_interval,
    interval_row,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_REPORT = ROOT / "robin_base_5041_55440_interval_report.json"
INTERVAL_BACKEND = ROOT / "test_ca_all_profile_interval_certificate.py"


@dataclass(frozen=True)
class MinimumRow:
    n: int
    sigma_n: int
    factorization: str
    gap_lower: Decimal
    gap_upper: Decimal


def sigma_divisor_sieve(limit: int) -> list[int]:
    """Compute sigma(n) exactly for every n <= limit by divisor summation."""

    values = [0] * (limit + 1)
    for divisor in range(1, limit + 1):
        for multiple in range(divisor, limit + 1, divisor):
            values[multiple] += divisor
    return values


def smallest_prime_factor_sieve(limit: int) -> list[int]:
    """Return an exact smallest-prime-factor table through limit."""

    spf = list(range(limit + 1))
    if limit >= 1:
        spf[1] = 1
    for prime in range(2, math.isqrt(limit) + 1):
        if spf[prime] != prime:
            continue
        for multiple in range(prime * prime, limit + 1, prime):
            if spf[multiple] == multiple:
                spf[multiple] = prime
    return spf


def factorization_and_sigma(n: int, spf: list[int]) -> tuple[str, int]:
    """Factor n with the SPF table and independently reconstruct sigma(n)."""

    remaining = n
    factors: list[str] = []
    sigma_n = 1
    while remaining > 1:
        prime = spf[remaining]
        exponent = 0
        prime_power = 1
        geometric_sum = 1
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
            prime_power *= prime
            geometric_sum += prime_power
        factors.append(f"{prime}^{exponent}")
        sigma_n *= geometric_sum
    return " ".join(factors), sigma_n


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Certify Robin's strict logarithmic gap for every integer from "
            "5041 through 55440 using exact sigma values and directed Decimal "
            "intervals."
        )
    )
    parser.add_argument("--start", type=int, default=5_041)
    parser.add_argument("--end", type=int, default=55_440)
    parser.add_argument("--precision", type=int, default=44)
    parser.add_argument("--gamma-terms", type=int, default=200_000)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.start < 3:
        raise ValueError("start must be at least 3")
    if args.end < args.start:
        raise ValueError("end must not be below start")

    started = time.perf_counter()
    arithmetic = DecimalIntervals(args.precision)
    gamma = harmonic_gamma_interval(arithmetic, args.gamma_terms)
    sigma_values = sigma_divisor_sieve(args.end)
    spf = smallest_prime_factor_sieve(args.end)

    checked = 0
    positive = 0
    sigma_crosschecks = 0
    minimum: MinimumRow | None = None
    maximum_width = Decimal(0)
    maximum_width_at: int | None = None
    failure_examples: list[dict[str, object]] = []

    for n in range(args.start, args.end + 1):
        factorization, sigma_factored = factorization_and_sigma(n, spf)
        sigma_sieved = sigma_values[n]
        if sigma_factored != sigma_sieved:
            raise AssertionError(
                f"exact sigma cross-check failed at n={n}: "
                f"{sigma_factored} != {sigma_sieved}"
            )
        sigma_crosschecks += 1

        log_n = arithmetic.ln_integer(n)
        log_log_n = arithmetic.ln(log_n)
        log_log_log_n = arithmetic.ln(log_log_n)
        sigma_over_n = arithmetic.rational(sigma_sieved, n)
        log_sigma_over_n = arithmetic.ln(sigma_over_n)
        gap = arithmetic.sub(
            arithmetic.add(gamma, log_log_log_n), log_sigma_over_n
        )

        checked += 1
        if gap[0] > 0:
            positive += 1
        elif len(failure_examples) < 25:
            failure_examples.append(
                {
                    "n": n,
                    "sigma_n": sigma_sieved,
                    "factorization": factorization,
                    "gap_interval": interval_row(gap),
                }
            )

        width = arithmetic.ceiling.subtract(gap[1], gap[0])
        if width > maximum_width:
            maximum_width = width
            maximum_width_at = n
        if minimum is None or gap[0] < minimum.gap_lower:
            minimum = MinimumRow(
                n=n,
                sigma_n=sigma_sieved,
                factorization=factorization,
                gap_lower=gap[0],
                gap_upper=gap[1],
            )

    assert minimum is not None
    interval_width_control = maximum_width < minimum.gap_lower
    overall_pass = (
        checked == args.end - args.start + 1
        and sigma_crosschecks == checked
        and positive == checked
        and not failure_examples
        and interval_width_control
    )

    script_path = Path(__file__).resolve()
    report = {
        "certificate_status": "PASS" if overall_pass else "FAIL",
        "claim_certified": (
            "For every integer n in the configured inclusive range, the "
            "strict logarithmic Robin gap gamma+log(log(log(n)))"
            "-log(sigma(n)/n) is positive."
        ),
        "configuration": {
            "start": args.start,
            "end": args.end,
            "precision_decimal_digits": args.precision,
            "gamma_harmonic_terms": args.gamma_terms,
        },
        "method": {
            "sigma_values": (
                "Exact Python integers from a divisor-sum sieve, independently "
                "cross-checked by SPF factorization and prime-power geometric sums."
            ),
            "strict_gap": (
                "For n>e^e, positivity of gamma+log(log(log(n)))"
                "-log(sigma(n)/n) is exactly equivalent to "
                "sigma(n)<exp(gamma)*n*log(log(n))."
            ),
            "directed_arithmetic": (
                "Every arithmetic operation is outward-rounded. Decimal ln is "
                "evaluated in the ROUND_HALF_EVEN context and padded by one ulp "
                "on each side. Euler's constant is enclosed from harmonic-number "
                "inequalities."
            ),
        },
        "coverage": {
            "integer_count": checked,
            "positive_gap_count": positive,
            "exact_sigma_crosscheck_count": sigma_crosschecks,
            "failure_count": checked - positive,
        },
        "euler_gamma_interval": interval_row(gamma),
        "minimum_gap": {
            "n": minimum.n,
            "sigma_n": minimum.sigma_n,
            "factorization": minimum.factorization,
            "gap_interval": {
                "lower": str(minimum.gap_lower),
                "upper": str(minimum.gap_upper),
            },
        },
        "interval_width_audit": {
            "maximum_gap_interval_width": str(maximum_width),
            "maximum_width_at_n": maximum_width_at,
            "maximum_width_is_below_uniform_positive_lower_bound": interval_width_control,
        },
        "failure_examples": failure_examples,
        "runtime_environment": {
            "python_version": sys.version,
            "python_implementation": platform.python_implementation(),
            "decimal_module_version": decimal.__version__,
            "libmpdec_version": decimal.__libmpdec_version__,
            "script_path": str(script_path),
            "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
            "interval_backend_path": str(INTERVAL_BACKEND.resolve()),
            "interval_backend_sha256": hashlib.sha256(
                INTERVAL_BACKEND.read_bytes()
            ).hexdigest(),
            "argv": sys.argv,
        },
        "timing_seconds": str(Decimal(str(time.perf_counter() - started))),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(
        f"{report['certificate_status']}: {checked:,} integers; "
        f"minimum lower gap {minimum.gap_lower} at n={minimum.n}; "
        f"report {args.report.resolve()}"
    )
    if not overall_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
