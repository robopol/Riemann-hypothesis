from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parent
DEFAULT_BLOCK_CSV = ROOT / "w_corrected_blocks_rows.csv"
DEFAULT_SUCCESSIVE_CSV = ROOT / "ca_successive_mvdc_blocks_rows.csv"
DEFAULT_JSON = ROOT / "w_chebyshev_target_report.json"
DEFAULT_CSV = ROOT / "w_chebyshev_target_rows.csv"


@dataclass
class ChebyshevTargetRow:
    y_prime: int
    x_prime: int
    f_gap: float
    s2_block: float
    corrected_w_block_exact: float
    chebyshev_integral: float
    c_max_exact: float
    c_max_m1: float
    sqrt_weight_integral: float
    average_k_from_integral: float
    required_k_for_c_190: float
    required_k_for_c_200: float
    required_k_for_c_205: float


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


def load_w_blocks(path: Path) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            rows.append(
                {
                    "y_prime": int(raw["y_prime"]),
                    "x_prime": int(raw["x_prime"]),
                    "w_block_exact": float(raw["w_block_exact"]),
                    "w_block_m1_upper": float(raw["w_block_m1_upper"]),
                }
            )
    return rows


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


def prime_block(primes: list[int], y: int, x: int) -> list[int]:
    y_index = nearest_prime_index_at_or_below(primes, y)
    x_index = nearest_prime_index_at_or_below(primes, x)
    return primes[y_index + 1 : x_index + 1]


def f_scale(x: int) -> float:
    return 1.0 / (math.sqrt(x) * math.log(x))


def s2_from_block(block: list[int]) -> float:
    return sum(-math.log1p(-1.0 / prime) - 1.0 / prime for prime in block)


def adaptive_simpson(
    fn: Callable[[float], float],
    left: float,
    right: float,
    eps: float = 1e-13,
    max_depth: int = 30,
) -> float:
    def simpson(a: float, b: float) -> float:
        c = 0.5 * (a + b)
        return (b - a) * (fn(a) + 4.0 * fn(c) + fn(b)) / 6.0

    def recurse(a: float, b: float, whole: float, depth: int) -> float:
        c = 0.5 * (a + b)
        left_part = simpson(a, c)
        right_part = simpson(c, b)
        total = left_part + right_part
        if depth <= 0 or abs(total - whole) <= 15.0 * eps:
            return total + (total - whole) / 15.0
        return recurse(a, c, left_part, depth - 1) + recurse(c, b, right_part, depth - 1)

    whole_interval = simpson(left, right)
    return recurse(left, right, whole_interval, max_depth)


def sqrt_weight_integral(y: int, x: int) -> float:
    # t=e^u transforms sqrt(t)*omega(t) dt into (u+1)/(sqrt(e^u)*u^2) du.
    left = math.log(y)
    right = math.log(x)

    def integrand(u: float) -> float:
        return (u + 1.0) * math.exp(-0.5 * u) / (u * u)

    return adaptive_simpson(integrand, left, right)


def required_k(c_value: float, f_gap: float, s2_block: float, sqrt_weight: float) -> float:
    return (c_value * f_gap + s2_block) / sqrt_weight


def evaluate(primes: list[int], block_row: dict[str, float | int]) -> ChebyshevTargetRow:
    y = int(block_row["y_prime"])
    x = int(block_row["x_prime"])
    block = prime_block(primes, y, x)
    s2 = s2_from_block(block)
    f_gap = f_scale(y) - f_scale(x)
    w_block = float(block_row["w_block_exact"])
    w_block_m1 = float(block_row["w_block_m1_upper"])
    cheb_integral = w_block - s2
    sqrt_weight = sqrt_weight_integral(y, x)
    return ChebyshevTargetRow(
        y_prime=y,
        x_prime=x,
        f_gap=f_gap,
        s2_block=s2,
        corrected_w_block_exact=w_block,
        chebyshev_integral=cheb_integral,
        c_max_exact=-w_block / f_gap,
        c_max_m1=-w_block_m1 / f_gap,
        sqrt_weight_integral=sqrt_weight,
        average_k_from_integral=-cheb_integral / sqrt_weight,
        required_k_for_c_190=required_k(1.90, f_gap, s2, sqrt_weight),
        required_k_for_c_200=required_k(2.00, f_gap, s2, sqrt_weight),
        required_k_for_c_205=required_k(2.05, f_gap, s2, sqrt_weight),
    )


def write_csv(path: Path, rows: list[ChebyshevTargetRow]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit the corrected W Chebyshev target on sampled CA blocks.")
    parser.add_argument("--block-csv", type=Path, default=DEFAULT_BLOCK_CSV)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    block_rows = load_w_blocks(args.block_csv)
    if not block_rows:
        raise SystemExit("No corrected W block rows found. Run test_w_corrected_blocks.py first.")
    max_prime = max(int(row["x_prime"]) for row in block_rows)
    primes = primes_up_to(max_prime)
    rows = [evaluate(primes, row) for row in block_rows]
    payload = {
        "notes": [
            "The exact target is integral_Y^x (theta(t)-t) omega(t) dt + S2(Y,x) <= -C_W(F(Y)-F(x)).",
            "The corrected W block is exactly the left side.",
            "average_k_from_integral is the weighted average in t-theta(t) ~= K sqrt(t) units.",
            "required_k_for_c_* is the weighted K sufficient for the selected C_W.",
        ],
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    write_csv(args.csv, rows)

    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print()
    print("Corrected W Chebyshev target:")
    for row in rows:
        print(
            f"{row.y_prime:>9}->{row.x_prime:<9} "
            f"Cmax={row.c_max_exact:.6f} "
            f"Kavg={row.average_k_from_integral:.6f} "
            f"Kreq(C=2)={row.required_k_for_c_200:.6f} "
            f"S2/Fgap={row.s2_block / row.f_gap:.6f}"
        )


if __name__ == "__main__":
    main()
