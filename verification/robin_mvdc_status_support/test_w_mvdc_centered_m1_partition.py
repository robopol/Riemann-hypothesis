from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from test_w_mvdc_certified_envelope import (
    DEFAULT_THETA_CSV,
    delta_rho_lower_bound,
    exact_block_values,
    f_scale,
    pi_lower,
    pi_upper,
    primes_up_to,
    prime_indices,
    s_bounds,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_JSON = ROOT / "w_mvdc_centered_m1_partition_report.json"
DEFAULT_CSV = ROOT / "w_mvdc_centered_m1_partition_rows.csv"


@dataclass
class CenteredM1PartitionRow:
    y_prime: int
    x_prime: int
    partitions: int
    endpoint_mode: str
    exact_m1: float
    certified_m1_upper: float
    exact_q: float
    q_upper_m1_plus_exact_tail: float
    delta_rho_lower: float
    w_upper_m1_plus_exact_tail: float
    exact_w_block: float
    exact_c_step: float
    c_step_certified: float | None
    certificate_margin: float


def read_theta_rows(path: Path) -> dict[int, dict[str, float]]:
    rows: dict[int, dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            x = int(raw["last_prime"])
            rows[x] = {key: float(value) for key, value in raw.items() if key != "last_prime"}
    return rows


def geometric_partition(y: int, x: int, count: int) -> list[tuple[int, int]]:
    if count < 1:
        raise ValueError("Partition count must be positive")
    log_y = math.log(y)
    log_x = math.log(x)
    raw = [int(round(math.exp(log_y + (log_x - log_y) * i / count))) for i in range(count + 1)]
    raw[0] = y
    raw[-1] = x

    intervals: list[tuple[int, int]] = []
    left = raw[0]
    for right in raw[1:]:
        if right <= left:
            continue
        intervals.append((left, right))
        left = right
    return intervals


def certified_m1_upper_by_partition(y: int, x: int, mu: float, partitions: int) -> float:
    constant_coeff = math.exp(-mu) - 1.0
    reciprocal_coeff = math.exp(-mu)
    total = 0.0
    for left, right in geometric_partition(y, x, partitions):
        s0_lower, _s0_upper = s_bounds(left, right, 0)
        _s1_lower, s1_upper = s_bounds(left, right, 1)
        total += constant_coeff * s0_lower + reciprocal_coeff * s1_upper
    return total


def exact_m1_from_block(primes: list[int], indices: dict[int, int], y: int, x: int, mu: float) -> float:
    y_index = indices[y]
    x_index = indices[x]
    base = math.exp(-mu) - 1.0
    coeff = math.exp(-mu)
    return sum(base + coeff / (prime - 1.0) for prime in primes[y_index + 1 : x_index + 1])


def exact_higher_taylor_tail(primes: list[int], indices: dict[int, int], y: int, x: int, mu: float) -> float:
    y_index = indices[y]
    x_index = indices[x]
    total = 0.0
    base = math.exp(-mu) - 1.0
    coeff = math.exp(-mu)
    for prime in primes[y_index + 1 : x_index + 1]:
        u = base + coeff / (prime - 1.0)
        total += math.log1p(u) - u
    return total


def evaluate_pair(
    primes: list[int],
    indices: dict[int, int],
    theta_rows: dict[int, dict[str, float]],
    y: int,
    x: int,
    partitions: int,
    endpoint_mode: str,
    theta_c: float,
    theta_power: int,
) -> CenteredM1PartitionRow:
    exact_nu, _main, q_value, _exact_s1, _exact_s2 = exact_block_values(primes, indices, y, x)
    main = math.log(math.log(x) / math.log(y))
    mu = main / exact_nu
    exact_m1 = exact_m1_from_block(primes, indices, y, x, mu)
    m1_upper = certified_m1_upper_by_partition(y, x, mu, partitions)
    higher_tail = exact_higher_taylor_tail(primes, indices, y, x, mu)
    q_upper = m1_upper + higher_tail
    delta_lower = delta_rho_lower_bound(y, x, theta_rows, endpoint_mode, theta_c, theta_power)
    w_upper = q_upper - delta_lower
    delta_exact = theta_rows[x]["endpoint_term"] - theta_rows[y]["endpoint_term"]
    w_exact = q_value - delta_exact
    f_gap = f_scale(y) - f_scale(x)
    c_cert = (-w_upper / f_gap) if w_upper < 0.0 else None
    c_exact = -w_exact / f_gap
    return CenteredM1PartitionRow(
        y_prime=y,
        x_prime=x,
        partitions=partitions,
        endpoint_mode=endpoint_mode,
        exact_m1=exact_m1,
        certified_m1_upper=m1_upper,
        exact_q=q_value,
        q_upper_m1_plus_exact_tail=q_upper,
        delta_rho_lower=delta_lower,
        w_upper_m1_plus_exact_tail=w_upper,
        exact_w_block=w_exact,
        exact_c_step=c_exact,
        c_step_certified=c_cert,
        certificate_margin=w_upper - w_exact,
    )


def write_csv(path: Path, rows: list[CenteredM1PartitionRow]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Test direct partitioned certificates for the centred MVDC first moment M1.")
    parser.add_argument("--theta-csv", type=Path, default=DEFAULT_THETA_CSV)
    parser.add_argument("--max-prime", type=int, default=120_000_000)
    parser.add_argument("--min-prime", type=int, default=3_329_267)
    parser.add_argument("--partitions", type=int, default=64)
    parser.add_argument("--endpoint-mode", choices=["exact", "theta-bound"], default="exact")
    parser.add_argument("--theta-c", type=float, default=0.006788)
    parser.add_argument("--theta-power", type=int, default=1)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    theta_rows = read_theta_rows(args.theta_csv)
    targets = [prime for prime in sorted(theta_rows) if args.min_prime <= prime <= args.max_prime]
    if len(targets) < 2:
        raise ValueError("Need at least two theta rows in range")

    primes = primes_up_to(max(args.max_prime, targets[-1]))
    indices = prime_indices(primes)
    rows = [
        evaluate_pair(
            primes=primes,
            indices=indices,
            theta_rows=theta_rows,
            y=y,
            x=x,
            partitions=args.partitions,
            endpoint_mode=args.endpoint_mode,
            theta_c=args.theta_c,
            theta_power=args.theta_power,
        )
        for y, x in zip(targets, targets[1:])
    ]

    payload = {
        "notes": [
            "This tests a direct certificate for the centred first moment M1.",
            "Each block is split geometrically. On every subblock the negative constant coefficient uses an S0 lower bound and the positive reciprocal coefficient uses an S1 upper bound.",
            "The higher Taylor contribution is kept exact in this diagnostic, so any remaining loss is the M1 certificate loss.",
            "The current implementation still fixes mu from the exact finite block count; a proof version needs a mu-interval certificate.",
        ],
        "args": {
            "theta_csv": str(args.theta_csv),
            "max_prime": args.max_prime,
            "min_prime": args.min_prime,
            "partitions": args.partitions,
            "endpoint_mode": args.endpoint_mode,
            "theta_c": args.theta_c,
            "theta_power": args.theta_power,
        },
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    write_csv(args.csv, rows)

    print(f"primes up to {max(args.max_prime, targets[-1])}: {len(primes)}")
    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print()
    print(f"Partitioned centred M1 certificate ({args.partitions} partitions, {args.endpoint_mode} endpoint):")
    for row in rows:
        cert = row.c_step_certified if row.c_step_certified is not None else float("nan")
        print(
            f"{row.y_prime:>9}->{row.x_prime:<9} "
            f"M1_exact={row.exact_m1:+.3e} M1_cert={row.certified_m1_upper:+.3e} "
            f"C_exact={row.exact_c_step:.6f} C_cert={cert:.6f} "
            f"margin={row.certificate_margin:+.3e}"
        )


if __name__ == "__main__":
    main()
