from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "ca_cumulative_mvdc_envelope_rows.csv"
DEFAULT_JSON = ROOT / "required_prime_harmonic_upper_report.json"
DEFAULT_CSV = ROOT / "required_prime_harmonic_upper_rows.csv"

MEISSEL_MERTENS_PRIME_CONSTANT = 0.261497212847642783755426838608695859


@dataclass
class RequiredPrimeHarmonicUpperRow:
    y_prime: int
    x_prime: int
    block_prime_count: int
    h_log: float
    mu: float
    actual_a_x: float
    required_a_upper: float
    actual_margin: float
    rosser_schoenfeld_finite_a_upper: float
    rosser_minus_required: float
    exact_moment_1: float
    m1_step_threshold: float


def primes_up_to(limit: int) -> list[int]:
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute the required upper envelope for A(x)=sum_{p<=x}1/p-log log x."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    with args.input.open(newline="", encoding="utf-8") as handle:
        input_rows = list(csv.DictReader(handle))
    if not input_rows:
        raise ValueError(f"No rows found in {args.input}")

    max_x = max(int(row["x_prime"]) for row in input_rows)
    primes = primes_up_to(max_x)

    prefix_inv_p = [0.0]
    prefix_c2 = [0.0]
    prefix_s_minus1 = [0.0]
    for prime in primes:
        prefix_inv_p.append(prefix_inv_p[-1] + 1.0 / prime)
        prefix_c2.append(prefix_c2[-1] + 1.0 / (prime * (prime - 1.0)))
        prefix_s_minus1.append(prefix_s_minus1[-1] + 1.0 / (prime - 1.0))

    def prime_index(value: int) -> int:
        return bisect.bisect_right(primes, value)

    def a_remainder(value: int) -> float:
        index = prime_index(value)
        return prefix_inv_p[index] - math.log(math.log(value))

    def block_c2(y: int, x: int) -> float:
        lo = prime_index(y)
        hi = prime_index(x)
        return prefix_c2[hi] - prefix_c2[lo]

    def block_s_minus1(y: int, x: int) -> float:
        lo = prime_index(y)
        hi = prime_index(x)
        return prefix_s_minus1[hi] - prefix_s_minus1[lo]

    def prime_count(y: int, x: int) -> int:
        return prime_index(x) - prime_index(y)

    def rosser_schoenfeld_finite_upper(x: int) -> float:
        return MEISSEL_MERTENS_PRIME_CONSTANT + 2.0 / (math.sqrt(x) * math.log(x))

    rows: list[RequiredPrimeHarmonicUpperRow] = []
    for raw in input_rows:
        y = int(raw["y_prime"])
        x = int(raw["x_prime"])
        nu = prime_count(y, x)
        h_log = math.log(math.log(x) / math.log(y))
        mu = h_log / nu
        exp_mu = math.exp(mu)
        c2 = block_c2(y, x)
        d_term = nu * (1.0 - (1.0 + mu) * math.exp(-mu))
        threshold = float(raw["m1_step_threshold"])

        # The step closes if A(x) is at most this value.
        required_a_upper = a_remainder(y) - c2 + exp_mu * (threshold + d_term)

        s_minus1 = block_s_minus1(y, x)
        exact_moment_1 = math.exp(-mu) * (nu + s_minus1) - nu
        actual_a_x = a_remainder(x)
        rs_upper = rosser_schoenfeld_finite_upper(x)

        rows.append(
            RequiredPrimeHarmonicUpperRow(
                y_prime=y,
                x_prime=x,
                block_prime_count=nu,
                h_log=h_log,
                mu=mu,
                actual_a_x=actual_a_x,
                required_a_upper=required_a_upper,
                actual_margin=required_a_upper - actual_a_x,
                rosser_schoenfeld_finite_a_upper=rs_upper,
                rosser_minus_required=rs_upper - required_a_upper,
                exact_moment_1=exact_moment_1,
                m1_step_threshold=threshold,
            )
        )

    payload = {
        "notes": [
            "This script computes the exact A(x) upper envelope required for the step M1 <= R_i-U_{i-1}.",
            "For a block (Y,x], the required condition is A(x) <= A(Y)-C2(Y,x)+exp(mu)*(threshold+D).",
            "The Rosser-Schoenfeld finite-range bound B+2/(sqrt(x) log x) is compared against that required upper envelope.",
            "Positive actual_margin means the exact block passes; positive rosser_minus_required means the RS finite bound is too loose for this step.",
        ],
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    with args.csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))

    print(f"read {args.input}")
    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print()
    print("Required prime-harmonic A(x) upper envelope:")
    for row in rows:
        print(
            f"  {row.y_prime:>9}->{row.x_prime:<9} "
            f"A={row.actual_a_x:.12e} "
            f"A_req={row.required_a_upper:.12e} "
            f"actual_margin={row.actual_margin:+.3e} "
            f"RS-A_req={row.rosser_minus_required:+.3e}"
        )


if __name__ == "__main__":
    main()
