from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_THETA_CSV = ROOT / "ca_theta_cancellation_rows.csv"
DEFAULT_JSON = ROOT / "w_mvdc_reduced_blocks_report.json"
DEFAULT_CSV = ROOT / "w_mvdc_reduced_blocks_rows.csv"


@dataclass
class CorrectedMVDCBlockRow:
    y_prime: int
    x_prime: int
    order: int
    block_prime_count: int
    mu: float
    max_abs_u: float
    q_exact: float
    q_upper_reduced: float
    q_upper_margin: float
    delta_rho: float
    w_block_exact: float
    w_block_upper_reduced: float
    w_upper_margin: float
    c_step_exact: float
    c_step_upper_reduced: float | None
    scaled_w_y: float
    scaled_w_x: float
    moment_1: float
    moment_2: float
    moment_3: float
    moment_4: float
    moment_5: float
    moment_6: float
    tail_remainder_bound: float
    s1_recip_p_minus_1: float
    s2_recip_p_minus_1: float


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


def prime_indices(primes: list[int]) -> dict[int, int]:
    return {prime: index for index, prime in enumerate(primes)}


def build_prefixes(primes: list[int], order: int) -> dict[int, list[float]]:
    prefixes: dict[int, list[float]] = {0: []}
    for s in range(1, order + 1):
        prefixes[s] = []

    totals = {s: 0.0 for s in range(0, order + 1)}
    for index, prime in enumerate(primes, start=1):
        totals[0] = float(index)
        prefixes[0].append(totals[0])
        for s in range(1, order + 1):
            totals[s] += 1.0 / ((prime - 1) ** s)
            prefixes[s].append(totals[s])
    return prefixes


def prefix_interval(prefix: list[float], left_exclusive_index: int, right_inclusive_index: int) -> float:
    if right_inclusive_index <= left_exclusive_index:
        return 0.0
    right = prefix[right_inclusive_index]
    left = prefix[left_exclusive_index]
    return right - left


def read_theta_rows(path: Path) -> dict[int, dict[str, float]]:
    rows: dict[int, dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            x = int(raw["last_prime"])
            rows[x] = {key: float(value) for key, value in raw.items() if key != "last_prime"}
    return rows


def f_scale(x: int) -> float:
    return 1.0 / (math.sqrt(x) * math.log(x))


def reduced_moment(moment_index: int, mu: float, reciprocal_sums: dict[int, float]) -> float:
    base = math.exp(-mu) - 1.0
    coeff = math.exp(-mu)
    total = 0.0
    for s in range(0, moment_index + 1):
        total += (
            math.comb(moment_index, s)
            * (base ** (moment_index - s))
            * (coeff**s)
            * reciprocal_sums[s]
        )
    return total


def partial_log_series(moments: dict[int, float], order: int) -> float:
    return sum(((-1.0) ** (r - 1)) * moments[r] / r for r in range(1, order + 1))


def block_remainder_bound(primes: list[int], y_index: int, x_index: int, mu: float, order: int) -> tuple[float, float]:
    values = []
    center = math.exp(mu)
    for prime in primes[y_index + 1 : x_index + 1]:
        values.append((1.0 / (1.0 - 1.0 / prime)) / center - 1.0)
    max_abs = max(abs(value) for value in values)
    total = 0.0
    for value in values:
        q = abs(value)
        if q >= 1.0:
            return math.inf, max_abs
        total += (q ** (order + 1)) / ((order + 1) * (1.0 - q))
    return total, max_abs


def q_exact(primes: list[int], y_index: int, x_index: int, main: float) -> float:
    return sum(-math.log1p(-1.0 / prime) for prime in primes[y_index + 1 : x_index + 1]) - main


def evaluate_pair(
    primes: list[int],
    indices: dict[int, int],
    prefixes: dict[int, list[float]],
    theta_rows: dict[int, dict[str, float]],
    y: int,
    x: int,
    order: int,
) -> CorrectedMVDCBlockRow:
    y_index = indices[y]
    x_index = indices[x]
    if x_index <= y_index:
        raise ValueError("Block endpoint order is invalid")

    reciprocal_sums: dict[int, float] = {}
    for s in range(0, order + 1):
        reciprocal_sums[s] = prefix_interval(prefixes[s], y_index, x_index)

    nu = int(round(reciprocal_sums[0]))
    main = math.log(math.log(x) / math.log(y))
    mu = main / nu
    moments = {r: reduced_moment(r, mu, reciprocal_sums) for r in range(1, order + 1)}
    remainder, max_abs_u = block_remainder_bound(primes, y_index, x_index, mu, order)
    q_upper = partial_log_series(moments, order) + remainder
    q_value = q_exact(primes, y_index, x_index, main)

    y_row = theta_rows[y]
    x_row = theta_rows[x]
    delta_rho = x_row["endpoint_term"] - y_row["endpoint_term"]
    w_exact = q_value - delta_rho
    w_upper = q_upper - delta_rho
    f_gap = f_scale(y) - f_scale(x)
    c_exact = -w_exact / f_gap
    c_upper = (-w_upper / f_gap) if w_upper < 0.0 and f_gap > 0.0 else None

    return CorrectedMVDCBlockRow(
        y_prime=y,
        x_prime=x,
        order=order,
        block_prime_count=nu,
        mu=mu,
        max_abs_u=max_abs_u,
        q_exact=q_value,
        q_upper_reduced=q_upper,
        q_upper_margin=q_upper - q_value,
        delta_rho=delta_rho,
        w_block_exact=w_exact,
        w_block_upper_reduced=w_upper,
        w_upper_margin=w_upper - w_exact,
        c_step_exact=c_exact,
        c_step_upper_reduced=c_upper,
        scaled_w_y=theta_rows[y]["scaled_modified_error_w"],
        scaled_w_x=theta_rows[x]["scaled_modified_error_w"],
        moment_1=moments.get(1, 0.0),
        moment_2=moments.get(2, 0.0),
        moment_3=moments.get(3, 0.0),
        moment_4=moments.get(4, 0.0),
        moment_5=moments.get(5, 0.0),
        moment_6=moments.get(6, 0.0),
        tail_remainder_bound=remainder,
        s1_recip_p_minus_1=reciprocal_sums.get(1, 0.0),
        s2_recip_p_minus_1=reciprocal_sums.get(2, 0.0),
    )


def write_csv(path: Path, rows: list[CorrectedMVDCBlockRow]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute corrected W-block MVDC upper bounds using reduced prime-block moments.")
    parser.add_argument("--theta-csv", type=Path, default=DEFAULT_THETA_CSV)
    parser.add_argument("--max-prime", type=int, default=120_000_000)
    parser.add_argument("--min-prime", type=int, default=3_329_267)
    parser.add_argument("--order", type=int, default=6)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.order < 1 or args.order > 6:
        raise ValueError("This report schema records moments only up to order 6")

    theta_rows = read_theta_rows(args.theta_csv)
    targets = [prime for prime in sorted(theta_rows) if prime >= args.min_prime and prime <= args.max_prime]
    if len(targets) < 2:
        raise ValueError("Need at least two theta rows in range")

    primes = primes_up_to(max(args.max_prime, targets[-1]))
    indices = prime_indices(primes)
    prefixes = build_prefixes(primes, args.order)

    rows = [
        evaluate_pair(primes, indices, prefixes, theta_rows, y, x, args.order)
        for y, x in zip(targets, targets[1:])
    ]

    payload = {
        "notes": [
            "Corrected block: W(x)-W(Y)=Q(Y,x)-(rho(x)-rho(Y)).",
            "Q is bounded by the MVDC Taylor polynomial using moments reduced to S_s=sum (p-1)^(-s).",
            "The remaining proof task is to replace the finite prime-block S_s and endpoint rho controls by certified analytic inequalities.",
        ],
        "args": {
            "theta_csv": str(args.theta_csv),
            "max_prime": args.max_prime,
            "min_prime": args.min_prime,
            "order": args.order,
        },
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    write_csv(args.csv, rows)

    print(f"primes up to {max(args.max_prime, targets[-1])}: {len(primes)}")
    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print()
    print("Corrected W MVDC reduced blocks:")
    for row in rows:
        c_upper = row.c_step_upper_reduced if row.c_step_upper_reduced is not None else float("nan")
        print(
            f"{row.y_prime:>9}->{row.x_prime:<9} "
            f"Wscale {row.scaled_w_y:.3f}->{row.scaled_w_x:.3f} "
            f"C_exact={row.c_step_exact:.6f} C_upper={c_upper:.6f} "
            f"Q_margin={row.q_upper_margin:+.2e} W_margin={row.w_upper_margin:+.2e}"
        )


if __name__ == "__main__":
    main()
