from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "ca_sminus_step_threshold_rows.csv"
DEFAULT_JSON = ROOT / "prime_harmonic_remainder_route_report.json"
DEFAULT_CSV = ROOT / "prime_harmonic_remainder_route_rows.csv"


@dataclass
class PrimeHarmonicRemainderRow:
    y_prime: int
    x_prime: int
    block_prime_count: int
    h_log: float
    mu: float
    harmonic_remainder_y: float
    harmonic_remainder_x: float
    harmonic_remainder_drop: float
    phi_from_remainders: float
    c2_correction: float
    convex_slack: float
    no_reserve_required_drop: float
    no_reserve_margin: float
    s_minus1: float
    s_minus1_threshold_no_reserve: float
    moment_1: float
    input_moment_1: float
    moment_1_difference: float


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
        description=(
            "Check the exact prime-harmonic remainder form of the S_{-1} block gate."
        )
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

    rows: list[PrimeHarmonicRemainderRow] = []
    for raw in input_rows:
        y = int(raw["y_prime"])
        x = int(raw["x_prime"])
        lo = bisect.bisect_right(primes, y)
        hi = bisect.bisect_right(primes, x)
        nu = hi - lo
        h_log = math.log(math.log(x) / math.log(y))
        mu = h_log / nu

        sum_inv_y = prefix_inv_p[lo]
        sum_inv_x = prefix_inv_p[hi]
        remainder_y = sum_inv_y - math.log(math.log(y))
        remainder_x = sum_inv_x - math.log(math.log(x))
        remainder_drop = remainder_y - remainder_x
        phi = remainder_x - remainder_y

        c2 = prefix_c2[hi] - prefix_c2[lo]
        s_minus1 = prefix_s_minus1[hi] - prefix_s_minus1[lo]
        convex_slack = nu * math.expm1(mu) - h_log
        no_reserve_required_drop = c2 - convex_slack
        no_reserve_margin = remainder_drop - no_reserve_required_drop
        threshold_no_reserve = nu * math.expm1(mu)
        moment_1 = math.exp(-mu) * (nu + s_minus1) - nu
        input_moment_1 = float(raw["moment_1"])

        rows.append(
            PrimeHarmonicRemainderRow(
                y_prime=y,
                x_prime=x,
                block_prime_count=nu,
                h_log=h_log,
                mu=mu,
                harmonic_remainder_y=remainder_y,
                harmonic_remainder_x=remainder_x,
                harmonic_remainder_drop=remainder_drop,
                phi_from_remainders=phi,
                c2_correction=c2,
                convex_slack=convex_slack,
                no_reserve_required_drop=no_reserve_required_drop,
                no_reserve_margin=no_reserve_margin,
                s_minus1=s_minus1,
                s_minus1_threshold_no_reserve=threshold_no_reserve,
                moment_1=moment_1,
                input_moment_1=input_moment_1,
                moment_1_difference=moment_1 - input_moment_1,
            )
        )

    payload = {
        "notes": [
            "Define A(x)=sum_{p<=x} 1/p - log log x.",
            "The exact identity is S_{-1}(Y,x)=H + A(x)-A(Y) + C2(Y,x).",
            "Without incoming reserve, the S_{-1} gate is equivalent to A(Y)-A(x) >= C2 - (nu*(exp(mu)-1)-H).",
            "Positive no_reserve_margin means the block already passes without using the cumulative reserve.",
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
    print("Prime-harmonic remainder gate:")
    for row in rows:
        print(
            f"  {row.y_prime:>9}->{row.x_prime:<9} "
            f"drop={row.harmonic_remainder_drop:.12e} "
            f"required={row.no_reserve_required_drop:.12e} "
            f"margin={row.no_reserve_margin:+.3e} "
            f"M1={row.moment_1:+.3e} "
            f"M1-diff={row.moment_1_difference:+.3e}"
        )


if __name__ == "__main__":
    main()
