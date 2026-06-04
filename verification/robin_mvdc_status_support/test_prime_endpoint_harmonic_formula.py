from __future__ import annotations

import argparse
import csv
import json
import math
import time
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_JSON = ROOT / "prime_endpoint_harmonic_formula_report.json"
DEFAULT_CSV = ROOT / "prime_endpoint_harmonic_formula_rows.csv"

MEISSEL_MERTENS_PRIME_CONSTANT = 0.261497212847642783755426838608695859


@dataclass
class PrimeEndpointRow:
    prime_index: int
    x_prime: int
    reciprocal_prime_sum: float
    loglog_plus_b: float
    formula_upper_c: float
    formula_margin: float
    formula_margin_percent_of_sum: float
    actual_remainder_minus_b: float
    actual_constant: float
    tested_constant: float
    constant_margin: float
    constant_margin_percent_of_actual: float | None
    passes_formula_upper: bool


def parse_float_list(value: str) -> list[float]:
    return [float(part.strip()) for part in value.split(",") if part.strip()]


def simple_primes_up_to(limit: int) -> list[int]:
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    root = math.isqrt(limit)
    for value in range(2, root + 1):
        if sieve[value]:
            start = value * value
            sieve[start : limit + 1 : value] = b"\x00" * (((limit - start) // value) + 1)
    return [value for value in range(2, limit + 1) if sieve[value]]


def segmented_primes(limit: int, segment_size: int):
    base_primes = simple_primes_up_to(math.isqrt(limit))
    if limit >= 2:
        yield 2

    for low in range(3, limit + 1, segment_size):
        high = min(limit, low + segment_size - 1)
        if low % 2 == 0:
            low += 1
        if low > high:
            continue

        size = ((high - low) // 2) + 1
        sieve = bytearray(b"\x01") * size
        for prime in base_primes[1:]:
            square = prime * prime
            if square > high:
                break
            start = max(square, ((low + prime - 1) // prime) * prime)
            if start % 2 == 0:
                start += prime
            if start > high:
                continue
            index = (start - low) // 2
            sieve[index::prime] = b"\x00" * (((size - 1 - index) // prime) + 1)

        for index, is_prime in enumerate(sieve):
            if is_prime:
                yield low + 2 * index


def evaluate_prime_endpoint(
    prime_index: int,
    x_prime: int,
    reciprocal_prime_sum: float,
    constant: float,
) -> PrimeEndpointRow:
    log_x = math.log(x_prime)
    scale = math.sqrt(x_prime) * log_x
    loglog_plus_b = math.log(log_x) + MEISSEL_MERTENS_PRIME_CONSTANT
    formula_upper_c = loglog_plus_b + constant / scale
    actual_remainder_minus_b = reciprocal_prime_sum - loglog_plus_b
    actual_constant = actual_remainder_minus_b * scale
    constant_margin = constant - actual_constant
    formula_margin = formula_upper_c - reciprocal_prime_sum
    constant_margin_percent = None
    if actual_constant != 0.0:
        constant_margin_percent = 100.0 * constant_margin / abs(actual_constant)
    return PrimeEndpointRow(
        prime_index=prime_index,
        x_prime=x_prime,
        reciprocal_prime_sum=reciprocal_prime_sum,
        loglog_plus_b=loglog_plus_b,
        formula_upper_c=formula_upper_c,
        formula_margin=formula_margin,
        formula_margin_percent_of_sum=100.0 * formula_margin / reciprocal_prime_sum,
        actual_remainder_minus_b=actual_remainder_minus_b,
        actual_constant=actual_constant,
        tested_constant=constant,
        constant_margin=constant_margin,
        constant_margin_percent_of_actual=constant_margin_percent,
        passes_formula_upper=formula_margin >= 0.0,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate the prime-endpoint formula sum_{j<=k}1/p_j <= "
            "log log p_k + B + C/(sqrt(p_k) log p_k)."
        )
    )
    parser.add_argument("--constant", type=float, default=1.5)
    parser.add_argument("--min-prime", type=int, default=2)
    parser.add_argument("--max-prime", type=int, default=1_000_000_000)
    parser.add_argument("--prime-index-step", type=int, default=5_000_000)
    parser.add_argument("--segment-size", type=int, default=8_000_000)
    parser.add_argument("--progress-every", type=int, default=100_000_000)
    parser.add_argument(
        "--scan-all-worst",
        action="store_true",
        help="Track the worst actual constant over every prime endpoint without writing every row.",
    )
    parser.add_argument(
        "--threshold-constants",
        default="",
        help="Comma-separated constants for eventual-threshold scanning. The main --constant is always included.",
    )
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.prime_index_step <= 0:
        raise ValueError("--prime-index-step must be positive.")

    rows: list[PrimeEndpointRow] = []
    reciprocal_prime_sum = 0.0
    prime_count = 0
    last_prime = 0
    worst_all: PrimeEndpointRow | None = None
    threshold_constants = sorted(set([args.constant] + parse_float_list(args.threshold_constants)))
    first_scanned: dict[float, PrimeEndpointRow | None] = {constant: None for constant in threshold_constants}
    last_failure: dict[float, PrimeEndpointRow | None] = {constant: None for constant in threshold_constants}
    first_after_last_failure: dict[float, PrimeEndpointRow | None] = {
        constant: None for constant in threshold_constants
    }
    start_time = time.time()
    last_progress = 0

    for prime in segmented_primes(args.max_prime, args.segment_size):
        reciprocal_prime_sum += 1.0 / prime
        prime_count += 1
        last_prime = prime

        if args.scan_all_worst and prime >= args.min_prime:
            candidate = evaluate_prime_endpoint(
                prime_count,
                prime,
                reciprocal_prime_sum,
                args.constant,
            )
            if worst_all is None or candidate.actual_constant > worst_all.actual_constant:
                worst_all = candidate
            for constant in threshold_constants:
                if first_scanned[constant] is None:
                    first_scanned[constant] = evaluate_prime_endpoint(
                        prime_count,
                        prime,
                        reciprocal_prime_sum,
                        constant,
                    )
                if candidate.actual_constant > constant:
                    last_failure[constant] = evaluate_prime_endpoint(
                        prime_count,
                        prime,
                        reciprocal_prime_sum,
                        constant,
                    )
                    first_after_last_failure[constant] = None
                elif first_after_last_failure[constant] is None:
                    first_after_last_failure[constant] = evaluate_prime_endpoint(
                        prime_count,
                        prime,
                        reciprocal_prime_sum,
                        constant,
                    )

        if prime_count % args.prime_index_step == 0:
            rows.append(
                evaluate_prime_endpoint(
                    prime_count,
                    prime,
                    reciprocal_prime_sum,
                    args.constant,
                )
            )

        if args.progress_every and prime >= last_progress + args.progress_every:
            elapsed = time.time() - start_time
            print(
                f"progress prime={prime:,} pi={prime_count:,} "
                f"elapsed={elapsed:.1f}s"
            )
            last_progress = prime

    if prime_count and (not rows or rows[-1].prime_index != prime_count):
        rows.append(
            evaluate_prime_endpoint(
                prime_count,
                last_prime,
                reciprocal_prime_sum,
                args.constant,
            )
        )

    failures = [row for row in rows if not row.passes_formula_upper]
    worst = min(rows, key=lambda row: row.constant_margin)
    summary_worst = worst_all if worst_all is not None else worst
    threshold_summary = {}
    for constant in threshold_constants:
        threshold_row = first_after_last_failure[constant]
        failure_row = last_failure[constant]
        if failure_row is None:
            threshold_row = first_scanned[constant]
        threshold_summary[str(constant)] = {
            "last_failure_prime_index": failure_row.prime_index if failure_row else None,
            "last_failure_x_prime": failure_row.x_prime if failure_row else None,
            "last_failure_actual_constant": failure_row.actual_constant if failure_row else None,
            "eventual_threshold_prime_index": threshold_row.prime_index if threshold_row else None,
            "eventual_threshold_x_prime": threshold_row.x_prime if threshold_row else None,
            "threshold_formula_margin_percent_of_sum": (
                threshold_row.formula_margin_percent_of_sum if threshold_row else None
            ),
            "threshold_constant_margin": threshold_row.constant_margin if threshold_row else None,
            "holds_through_max_prime_after_threshold": threshold_row is not None,
        }
    payload = {
        "notes": [
            "Numerical audit only; this is not an analytic proof.",
            "The endpoint is x=p_k, not a real fixed grid point.",
            "B is the Meissel-Mertens prime constant.",
            "The tested formula is sum_{j<=k}1/p_j <= log log p_k + B + C/(sqrt(p_k) log p_k).",
            "actual_constant = (sum_{j<=k}1/p_j - log log p_k - B)*sqrt(p_k)*log(p_k).",
            "The formula with constant C passes exactly when actual_constant <= C.",
            "formula_margin_percent_of_sum = 100*(formula upper - exact sum)/exact sum.",
            "For each threshold constant, eventual_threshold_x_prime is the first prime endpoint after the last failure within the scanned range.",
        ],
        "args": {
            "constant": args.constant,
            "min_prime": args.min_prime,
            "max_prime": args.max_prime,
            "prime_index_step": args.prime_index_step,
            "segment_size": args.segment_size,
            "threshold_constants": threshold_constants,
        },
        "summary": {
            "tested_prime_endpoints": len(rows),
            "failures": len(failures),
            "worst_prime_index": worst.prime_index,
            "worst_x_prime": worst.x_prime,
            "worst_actual_constant": worst.actual_constant,
            "worst_constant_margin": worst.constant_margin,
            "scan_all_worst_enabled": args.scan_all_worst,
            "all_prime_endpoints_scanned": prime_count if args.scan_all_worst else None,
            "all_worst_min_prime": args.min_prime if args.scan_all_worst else None,
            "all_worst_prime_index": summary_worst.prime_index,
            "all_worst_x_prime": summary_worst.x_prime,
            "all_worst_actual_constant": summary_worst.actual_constant,
            "all_worst_constant_margin": summary_worst.constant_margin,
            "thresholds": threshold_summary,
        },
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    with args.csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))

    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print()
    print(
        f"Prime-endpoint formula for C={args.constant:g}: "
        f"{len(rows) - len(failures)}/{len(rows)} endpoints pass"
    )
    print(
        f"worst k={worst.prime_index:,} x=p_k={worst.x_prime:,} "
        f"C_actual={worst.actual_constant:.6f} "
        f"C_margin={worst.constant_margin:+.6f}"
    )
    if worst_all is not None:
        print(
            f"all-prime scan worst k={worst_all.prime_index:,} "
            f"x=p_k={worst_all.x_prime:,} "
            f"C_actual={worst_all.actual_constant:.6f} "
            f"C_margin={worst_all.constant_margin:+.6f}"
        )
        print("eventual thresholds:")
        for constant in threshold_constants:
            threshold_row = first_after_last_failure[constant]
            failure_row = last_failure[constant]
            if failure_row is None:
                threshold_row = first_scanned[constant]
            if threshold_row is None:
                print(f"  C={constant:g}: no threshold inside scanned range")
                continue
            last_failure_text = (
                "none"
                if failure_row is None
                else f"k={failure_row.prime_index:,}, p_k={failure_row.x_prime:,}, C_actual={failure_row.actual_constant:.6f}"
            )
            print(
                f"  C={constant:g}: last_failure={last_failure_text}; "
                f"holds from k={threshold_row.prime_index:,}, p_k={threshold_row.x_prime:,}; "
                f"margin={threshold_row.formula_margin_percent_of_sum:.9f}% of sum"
            )
    for row in rows:
        print(
            f"  k={row.prime_index:>10,} p_k={row.x_prime:>12,} "
            f"sum1p={row.reciprocal_prime_sum:.12f} "
            f"formula_C={row.formula_upper_c:.12f} "
            f"margin_pct={row.formula_margin_percent_of_sum:+.9f}% "
            f"C_actual={row.actual_constant:+.6f} "
            f"C_margin={row.constant_margin:+.6f} "
            f"pass={row.passes_formula_upper}"
        )


if __name__ == "__main__":
    main()
