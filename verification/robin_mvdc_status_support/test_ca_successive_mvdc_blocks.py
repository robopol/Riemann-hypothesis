from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_CA_CSV = ROOT / "mvdc_ca_beta_target_high_alpha_ca_rows.csv"
DEFAULT_JSON = ROOT / "ca_successive_mvdc_blocks_report.json"
DEFAULT_CSV = ROOT / "ca_successive_mvdc_blocks_rows.csv"


@dataclass
class CASuccessiveBlockRow:
    y_prime: int
    x_prime: int
    block_prime_count: int
    q_exact: float
    moment_1: float
    moment_2: float
    q_minus_m1: float
    c_step_q: float | None
    c_step_m1: float | None
    margin_q_c_130: float | None
    margin_q_c_135: float | None
    margin_q_c_140: float | None
    margin_m1_c_140: float | None


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


def load_ca_targets(path: Path, max_prime: int) -> list[int]:
    targets: list[int] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            prime = int(row["last_prime"])
            if prime <= max_prime:
                targets.append(prime)
    return sorted(set(targets))


def nearest_prime_index_at_or_below(primes: list[int], target: int) -> int:
    lo = 0
    hi = len(primes)
    while lo < hi:
        mid = (lo + hi) // 2
        if primes[mid] <= target:
            lo = mid + 1
        else:
            hi = mid
    if lo == 0:
        raise ValueError(f"No prime <= {target}")
    return lo - 1


def evaluate_pair(primes: list[int], y: int, x: int) -> CASuccessiveBlockRow:
    y_index = nearest_prime_index_at_or_below(primes, y)
    x_index = nearest_prime_index_at_or_below(primes, x)
    block = primes[y_index + 1 : x_index + 1]
    if not block:
        raise ValueError("Empty CA block")

    main = math.log(math.log(x) / math.log(y))
    q_exact = sum(-math.log1p(-1.0 / prime) for prime in block) - main
    center_h = math.exp(main / len(block))
    u_values = [(1.0 / (1.0 - 1.0 / prime)) / center_h - 1.0 for prime in block]
    m1 = sum(u_values)
    m2 = sum(value * value for value in u_values)
    f_gap = 1.0 / (math.sqrt(y) * math.log(y)) - 1.0 / (math.sqrt(x) * math.log(x))

    def c_step(value: float) -> float | None:
        return (-value / f_gap) if value < 0.0 and f_gap > 0.0 else None

    c_q = c_step(q_exact)
    c_m1 = c_step(m1)

    def margin(value: float | None, candidate: float) -> float | None:
        if value is None:
            return None
        return value - candidate

    return CASuccessiveBlockRow(
        y_prime=y,
        x_prime=x,
        block_prime_count=len(block),
        q_exact=q_exact,
        moment_1=m1,
        moment_2=m2,
        q_minus_m1=q_exact - m1,
        c_step_q=c_q,
        c_step_m1=c_m1,
        margin_q_c_130=margin(c_q, 1.30),
        margin_q_c_135=margin(c_q, 1.35),
        margin_q_c_140=margin(c_q, 1.40),
        margin_m1_c_140=margin(c_m1, 1.40),
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Test MVDC induction blocks between successive sampled CA support primes.")
    parser.add_argument("--max-prime", type=int, default=100_000_000)
    parser.add_argument("--ca-csv", type=Path, default=DEFAULT_CA_CSV)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    targets = load_ca_targets(args.ca_csv, args.max_prime)
    primes = primes_up_to(max(args.max_prime, max(targets)))

    rows = [evaluate_pair(primes, y, x) for y, x in zip(targets, targets[1:])]
    payload = {
        "notes": [
            "These are direct MVDC blocks between successive sampled CA support primes.",
            "If these steps pass for C, an induction can stay on the CA support instead of requiring an all-prime envelope at intermediate cutoffs.",
            "Q<=M1 is always available from log(1+u)<=u; the columns also report the M1-only step.",
        ],
        "ca_targets": targets,
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    with args.csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))

    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print()
    print("Successive CA MVDC blocks:")
    for row in rows:
        print(
            f"  {row.y_prime:>9}->{row.x_prime:<9} "
            f"Q={row.q_exact:+.6e} M1={row.moment_1:+.6e} "
            f"Cq={row.c_step_q if row.c_step_q is not None else float('nan'):.6f} "
            f"Cm1={row.c_step_m1 if row.c_step_m1 is not None else float('nan'):.6f} "
            f"margin C=1.40={row.margin_q_c_140 if row.margin_q_c_140 is not None else float('nan'):+.6f}"
        )


if __name__ == "__main__":
    main()
