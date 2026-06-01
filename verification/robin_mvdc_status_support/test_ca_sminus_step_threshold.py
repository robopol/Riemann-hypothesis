from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_CUMULATIVE_CSV = ROOT / "ca_cumulative_mvdc_envelope_rows.csv"
DEFAULT_JSON = ROOT / "ca_sminus_step_threshold_report.json"
DEFAULT_CSV = ROOT / "ca_sminus_step_threshold_rows.csv"


@dataclass
class CASMinusStepThresholdRow:
    y_prime: int
    x_prime: int
    block_prime_count: int
    h_log: float
    mu: float
    actual_s_minus1: float
    threshold_s_minus1: float
    margin_s_minus1: float
    moment_1: float
    m1_step_threshold: float
    m1_step_margin: float
    s_margin_from_m1_margin: float
    margin_difference: float


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
        description="Check the exact S_{-1} step threshold equivalent to the cumulative CA M1 ledger."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_CUMULATIVE_CSV)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    with args.input.open(newline="", encoding="utf-8") as handle:
        input_rows = list(csv.DictReader(handle))
    if not input_rows:
        raise ValueError(f"No rows in {args.input}")

    max_x = max(int(row["x_prime"]) for row in input_rows)
    primes = primes_up_to(max_x)

    rows: list[CASMinusStepThresholdRow] = []
    for raw in input_rows:
        y = int(raw["y_prime"])
        x = int(raw["x_prime"])
        block = [prime for prime in primes if y < prime <= x]
        nu = len(block)
        h_log = math.log(math.log(x) / math.log(y))
        mu = h_log / nu
        actual_s = sum(1.0 / (prime - 1.0) for prime in block)

        reserve_minus_previous_upper = float(raw["m1_step_threshold"])
        threshold = nu * math.expm1(mu) + math.exp(mu) * reserve_minus_previous_upper
        s_margin = threshold - actual_s

        m1 = float(raw["moment_1"])
        m1_margin = float(raw["m1_step_margin"])
        s_margin_from_m1_margin = math.exp(mu) * m1_margin

        rows.append(
            CASMinusStepThresholdRow(
                y_prime=y,
                x_prime=x,
                block_prime_count=nu,
                h_log=h_log,
                mu=mu,
                actual_s_minus1=actual_s,
                threshold_s_minus1=threshold,
                margin_s_minus1=s_margin,
                moment_1=m1,
                m1_step_threshold=reserve_minus_previous_upper,
                m1_step_margin=m1_margin,
                s_margin_from_m1_margin=s_margin_from_m1_margin,
                margin_difference=s_margin - s_margin_from_m1_margin,
            )
        )

    payload = {
        "notes": [
            "The checked threshold is S_{-1}(Y,x) <= nu*(exp(mu)-1) + exp(mu)*(R_i-U_{i-1}).",
            "It is algebraically equivalent to M1(Y,x) <= R_i-U_{i-1}.",
            "The final column verifies that the S-margin equals exp(mu) times the M1-margin up to floating-point noise.",
        ],
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    with args.csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))

    print(f"read {args.input}")
    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print()
    print("CA S_{-1} step threshold:")
    for row in rows:
        print(
            f"  {row.y_prime:>9}->{row.x_prime:<9} "
            f"S={row.actual_s_minus1:.12e} "
            f"T={row.threshold_s_minus1:.12e} "
            f"S-margin={row.margin_s_minus1:+.3e} "
            f"M1-margin={row.m1_step_margin:+.3e} "
            f"diff={row.margin_difference:+.3e}"
        )


if __name__ == "__main__":
    main()
