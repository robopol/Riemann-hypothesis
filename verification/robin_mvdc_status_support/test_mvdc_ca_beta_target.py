from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_JSON = ROOT / "mvdc_ca_beta_target_report.json"
DEFAULT_CA_CSV = ROOT / "mvdc_ca_beta_target_ca_rows.csv"
DEFAULT_BLOCK_CSV = ROOT / "mvdc_ca_beta_target_block_rows.csv"

EULER_GAMMA = 0.577215664901532860606512090082402431


@dataclass
class CABetaTargetRow:
    epsilon: float
    last_prime: int
    prime_count: int
    beta_error_e: float
    true_deficit_a: float
    bridge_b: float
    total_a_plus_b: float
    margin_vs_e: float
    c_beta_actual: float
    c_reserve_actual: float
    c_margin_actual: float
    first_unit_prime: int
    first_unit_over_sqrt_last: float
    last_nonunit_prime: int | None
    max_exponent: int
    log_n_minus_last_prime: float


@dataclass
class MVDCBlockRow:
    last_prime: int
    cutoff_y: int
    alpha: float
    block_prime_count: int
    beta_error_e_x: float
    beta_error_e_y: float
    q_exact: float
    partial_order: int
    partial_value: float
    remainder_bound: float
    q_upper: float
    q_upper_margin: float
    max_abs_u: float
    scaled_q: float
    scaled_q_upper: float
    induction_f_gap: float
    induction_c_step_max: float | None
    induction_margin_c_130: float | None
    induction_margin_c_135: float | None
    induction_margin_c_140: float | None
    induction_margin_c_145: float | None


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


def parse_float_list(value: str) -> list[float]:
    return [float(part.strip()) for part in value.split(",") if part.strip()]


def build_log_beta_by_prime(primes: list[int]) -> dict[int, float]:
    total = 0.0
    values: dict[int, float] = {}
    for prime in primes:
        total += -math.log1p(-1.0 / prime)
        values[prime] = total
    return values


def beta_error(log_beta: float, x: int) -> float:
    return log_beta - EULER_GAMMA - math.log(math.log(x))


def ca_exponent(prime: int, epsilon: float) -> int:
    log_p = math.log(prime)
    numerator_log = (1.0 + epsilon) * log_p + math.log1p(-math.exp(-(1.0 + epsilon) * log_p))
    denominator_log = math.log(math.expm1(epsilon * log_p))
    exponent = math.floor((numerator_log - denominator_log) / log_p) - 1
    return max(0, exponent)


def ca_profile(all_primes: list[int], epsilon: float) -> dict[int, int]:
    profile: dict[int, int] = {}
    for prime in all_primes:
        exponent = ca_exponent(prime, epsilon)
        if exponent > 0:
            profile[prime] = exponent
        elif profile:
            break
    return profile


def evaluate_ca_profile(
    epsilon: float,
    profile: dict[int, int],
    log_beta_by_prime: dict[int, float],
) -> CABetaTargetRow:
    primes = sorted(profile)
    last_prime = primes[-1]
    log_n = 0.0
    deficit = 0.0
    max_exponent = 0
    first_unit = None
    last_nonunit = None

    for prime in primes:
        exponent = profile[prime]
        max_exponent = max(max_exponent, exponent)
        log_n += exponent * math.log(prime)
        deficit += -math.log1p(-(prime ** (-(exponent + 1))))
        if exponent == 1 and first_unit is None:
            first_unit = prime
        if exponent > 1:
            last_nonunit = prime

    if first_unit is None:
        raise ValueError("CA profile has no unit tail in the generated range.")

    bridge = math.log(math.log(log_n) / math.log(last_prime)) if log_n > last_prime else 0.0
    error = beta_error(log_beta_by_prime[last_prime], last_prime)
    scale = math.sqrt(last_prime) * math.log(last_prime)
    total = deficit + bridge

    return CABetaTargetRow(
        epsilon=epsilon,
        last_prime=last_prime,
        prime_count=len(primes),
        beta_error_e=error,
        true_deficit_a=deficit,
        bridge_b=bridge,
        total_a_plus_b=total,
        margin_vs_e=total - error,
        c_beta_actual=error * scale,
        c_reserve_actual=total * scale,
        c_margin_actual=(total - error) * scale,
        first_unit_prime=first_unit,
        first_unit_over_sqrt_last=first_unit / math.sqrt(last_prime),
        last_nonunit_prime=last_nonunit,
        max_exponent=max_exponent,
        log_n_minus_last_prime=log_n - last_prime,
    )


def nearest_prime_at_or_below(primes: list[int], target: int) -> int:
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
    return primes[lo - 1]


def partial_log_series(moments: dict[int, float], order: int) -> float:
    return sum(((-1.0) ** (r - 1)) * moments[r] / r for r in range(1, order + 1))


def taylor_log_remainder_bound(u_values: list[float], order: int) -> float:
    total = 0.0
    for value in u_values:
        q = abs(value)
        if q >= 1.0:
            return math.inf
        total += (q ** (order + 1)) / ((order + 1) * (1.0 - q))
    return total


def evaluate_mvdc_block(
    primes: list[int],
    log_beta_by_prime: dict[int, float],
    x: int,
    alpha: float,
    order: int,
) -> MVDCBlockRow | None:
    y = nearest_prime_at_or_below(primes, int(x**alpha))
    if y >= x:
        return None
    block = [prime for prime in primes if y < prime <= x]
    if not block:
        return None

    e_x = beta_error(log_beta_by_prime[x], x)
    e_y = beta_error(log_beta_by_prime[y], y)
    main = math.log(math.log(x) / math.log(y))
    q_exact = sum(-math.log1p(-1.0 / prime) for prime in block) - main

    center_h = math.exp(main / len(block))
    u_values = [(1.0 / (1.0 - 1.0 / prime)) / center_h - 1.0 for prime in block]
    moments = {r: sum(value**r for value in u_values) for r in range(1, order + 1)}
    partial = partial_log_series(moments, order)
    remainder = taylor_log_remainder_bound(u_values, order)
    q_upper = partial + remainder
    scale = math.sqrt(x) * math.log(x)
    f_x = 1.0 / (math.sqrt(x) * math.log(x))
    f_y = 1.0 / (math.sqrt(y) * math.log(y))
    f_gap = f_y - f_x
    c_step_max = (-q_upper / f_gap) if q_upper < 0.0 and f_gap > 0.0 else None

    def induction_margin(candidate: float) -> float | None:
        if c_step_max is None:
            return None
        return c_step_max - candidate

    return MVDCBlockRow(
        last_prime=x,
        cutoff_y=y,
        alpha=alpha,
        block_prime_count=len(block),
        beta_error_e_x=e_x,
        beta_error_e_y=e_y,
        q_exact=q_exact,
        partial_order=order,
        partial_value=partial,
        remainder_bound=remainder,
        q_upper=q_upper,
        q_upper_margin=q_upper - q_exact,
        max_abs_u=max(abs(value) for value in u_values),
        scaled_q=q_exact * scale,
        scaled_q_upper=q_upper * scale,
        induction_f_gap=f_gap,
        induction_c_step_max=c_step_max,
        induction_margin_c_130=induction_margin(1.30),
        induction_margin_c_135=induction_margin(1.35),
        induction_margin_c_140=induction_margin(1.40),
        induction_margin_c_145=induction_margin(1.45),
    )


def scan_all_prime_beta_constant(primes: list[int], log_beta_by_prime: dict[int, float]) -> dict[str, float | int]:
    max_c = -math.inf
    max_prime = 0
    max_e = -math.inf
    max_e_prime = 0
    for prime in primes:
        if prime < 3:
            continue
        error = beta_error(log_beta_by_prime[prime], prime)
        scale = math.sqrt(prime) * math.log(prime)
        c_value = error * scale
        if c_value > max_c:
            max_c = c_value
            max_prime = prime
        if error > max_e:
            max_e = error
            max_e_prime = prime
    return {
        "max_c_beta_all_primes": max_c,
        "max_c_beta_prime": max_prime,
        "max_e_all_primes": max_e,
        "max_e_prime": max_e_prime,
    }


def scan_beta_constants_by_threshold(
    primes: list[int],
    log_beta_by_prime: dict[int, float],
    thresholds: list[int],
) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for threshold in thresholds:
        max_c = -math.inf
        max_prime = 0
        max_e = -math.inf
        max_e_prime = 0
        count = 0
        for prime in primes:
            if prime < max(3, threshold):
                continue
            count += 1
            error = beta_error(log_beta_by_prime[prime], prime)
            scale = math.sqrt(prime) * math.log(prime)
            c_value = error * scale
            if c_value > max_c:
                max_c = c_value
                max_prime = prime
            if error > max_e:
                max_e = error
                max_e_prime = prime
        rows.append(
            {
                "threshold": threshold,
                "prime_count": count,
                "max_c_beta": max_c,
                "max_c_beta_prime": max_prime,
                "max_e": max_e,
                "max_e_prime": max_e_prime,
            }
        )
    return rows


def write_csv(path: Path, rows: list[object]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Measure the MVDC beta-error target against the exact CA-profile ledger."
    )
    parser.add_argument("--max-prime", type=int, default=50_000_000)
    parser.add_argument(
        "--epsilons",
        default=(
            "0.001,0.0005,0.0002,0.0001,0.00005,0.00002,0.00001,"
            "0.000005,0.000002,0.000001,0.0000005,0.0000002,"
            "0.0000001,0.00000005,0.00000002,0.00000001,"
            "0.000000005,0.000000002,0.000000001"
        ),
    )
    parser.add_argument("--alphas", default="0.5,0.6666667,0.75,0.9,0.95")
    parser.add_argument("--order", type=int, default=8)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--ca-csv", type=Path, default=DEFAULT_CA_CSV)
    parser.add_argument("--block-csv", type=Path, default=DEFAULT_BLOCK_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    epsilons = parse_float_list(args.epsilons)
    alphas = parse_float_list(args.alphas)

    primes = primes_up_to(args.max_prime)
    log_beta_by_prime = build_log_beta_by_prime(primes)

    ca_rows: list[CABetaTargetRow] = []
    for epsilon in epsilons:
        profile = ca_profile(primes, epsilon)
        if profile:
            ca_rows.append(evaluate_ca_profile(epsilon, profile, log_beta_by_prime))

    block_rows: list[MVDCBlockRow] = []
    ca_last_primes = sorted({row.last_prime for row in ca_rows})
    for last_prime in ca_last_primes:
        for alpha in alphas:
            block_row = evaluate_mvdc_block(primes, log_beta_by_prime, last_prime, alpha, args.order)
            if block_row is not None:
                block_rows.append(block_row)

    ca_summary = {
        "min_c_reserve_actual": min(row.c_reserve_actual for row in ca_rows),
        "min_c_reserve_epsilon": min(ca_rows, key=lambda row: row.c_reserve_actual).epsilon,
        "min_c_reserve_prime": min(ca_rows, key=lambda row: row.c_reserve_actual).last_prime,
        "max_c_beta_actual_on_ca": max(row.c_beta_actual for row in ca_rows),
        "max_c_beta_epsilon": max(ca_rows, key=lambda row: row.c_beta_actual).epsilon,
        "max_c_beta_prime": max(ca_rows, key=lambda row: row.c_beta_actual).last_prime,
        "min_c_margin_actual": min(row.c_margin_actual for row in ca_rows),
        "min_margin_prime": min(ca_rows, key=lambda row: row.c_margin_actual).last_prime,
    }
    all_prime_scan = scan_all_prime_beta_constant(primes, log_beta_by_prime)
    threshold_scan = scan_beta_constants_by_threshold(
        primes,
        log_beta_by_prime,
        [1_000, 10_000, 100_000, 1_000_000, 2_000_000, 5_000_000, 10_000_000, 50_000_000],
    )
    best_induction_by_prime = {}
    for row in block_rows:
        current = best_induction_by_prime.get(row.last_prime)
        if row.induction_c_step_max is None:
            continue
        if current is None or row.induction_c_step_max > current["best_c_step_max"]:
            best_induction_by_prime[row.last_prime] = {
                "best_alpha": row.alpha,
                "best_cutoff_y": row.cutoff_y,
                "best_c_step_max": row.induction_c_step_max,
                "margin_c_130": row.induction_c_step_max - 1.30,
                "margin_c_135": row.induction_c_step_max - 1.35,
                "margin_c_140": row.induction_c_step_max - 1.40,
                "margin_c_145": row.induction_c_step_max - 1.45,
            }

    induction_summary = {
        "description": (
            "For F_C(x)=C/(sqrt(x)log(x)), one MVDC induction step from Y to x "
            "is certified when q_upper <= -C*(F_1(Y)-F_1(x)). "
            "Thus induction_c_step_max=-q_upper/(F_1(Y)-F_1(x))."
        ),
        "best_by_prime": best_induction_by_prime,
    }
    for candidate in (1.30, 1.35, 1.40, 1.45):
        passed = [
            row
            for row in block_rows
            if row.induction_c_step_max is not None and row.induction_c_step_max >= candidate
        ]
        failed = [
            row
            for row in block_rows
            if row.induction_c_step_max is None or row.induction_c_step_max < candidate
        ]
        induction_summary[f"C_{candidate:.2f}"] = {
            "passed_blocks": len(passed),
            "failed_blocks": len(failed),
            "best_alpha_passed_primes": sorted(
                prime
                for prime, best in best_induction_by_prime.items()
                if best["best_c_step_max"] >= candidate
            ),
            "best_alpha_failed_primes": sorted(
                prime
                for prime, best in best_induction_by_prime.items()
                if best["best_c_step_max"] < candidate
            ),
        }

    payload = {
        "notes": [
            "E(x)=log(beta(x))-gamma-log(log(x)).",
            "The target envelope tested here is E(x)<=C/(sqrt(x) log x).",
            "A CA ledger closes numerically when C is below c_reserve_actual=(A_CA+B_log_CA)*sqrt(x)*log(x).",
            "MVDC block rows give the exact renormalised block Q(Y,x) and a finite Taylor upper bound using exact prime sums.",
            "The induction columns test whether one block can propagate E(t)<=C/(sqrt(t)log(t)) from Y to x.",
            "This report is numerical. A rigorous proof still needs analytic upper bounds for the prime-block moments or another proof of the same C envelope.",
        ],
        "args": {
            "max_prime": args.max_prime,
            "epsilons": epsilons,
            "alphas": alphas,
            "order": args.order,
        },
        "ca_summary": ca_summary,
        "all_prime_scan": all_prime_scan,
        "threshold_scan": threshold_scan,
        "induction_summary": induction_summary,
        "candidate_constants": [
            {
                "C": c_value,
                "ca_margin_to_min_reserve": ca_summary["min_c_reserve_actual"] - c_value,
            }
            for c_value in (1.30, 1.35, 1.40, 1.45)
        ],
        "ca_rows": [asdict(row) for row in ca_rows],
        "block_rows": [asdict(row) for row in block_rows],
    }

    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    write_csv(args.ca_csv, ca_rows)
    write_csv(args.block_csv, block_rows)

    print(f"primes up to {args.max_prime}: {len(primes)}")
    print(f"wrote {args.json}")
    print(f"wrote {args.ca_csv}")
    print(f"wrote {args.block_csv}")
    print()
    print("CA beta target summary:")
    print(
        f"  min C_reserve on CA rows = {ca_summary['min_c_reserve_actual']:.6f} "
        f"at eps={ca_summary['min_c_reserve_epsilon']:g}, pk={ca_summary['min_c_reserve_prime']}"
    )
    print(
        f"  max C_beta on CA rows    = {ca_summary['max_c_beta_actual_on_ca']:.6f} "
        f"at eps={ca_summary['max_c_beta_epsilon']:g}, pk={ca_summary['max_c_beta_prime']}"
    )
    print(
        f"  max C_beta on all primes = {all_prime_scan['max_c_beta_all_primes']:.6f} "
        f"at pk={all_prime_scan['max_c_beta_prime']}"
    )
    for candidate in payload["candidate_constants"]:
        print(
            f"  C={candidate['C']:.2f}: margin to min CA reserve "
            f"{candidate['ca_margin_to_min_reserve']:+.6f}"
        )
    print()
    print("Suffix scans for C_beta=E(x)*sqrt(x)*log(x):")
    for row in threshold_scan:
        print(
            f"  p>={row['threshold']:>8}: max C_beta={row['max_c_beta']:.6f} "
            f"at pk={row['max_c_beta_prime']}"
        )
    print()
    print("Best MVDC induction step by CA prime:")
    for prime, best in sorted(best_induction_by_prime.items()):
        print(
            f"  pk={prime:>9} best alpha={best['best_alpha']:.3f} "
            f"y={best['best_cutoff_y']:>9} C_step_max={best['best_c_step_max']:.6f} "
            f"margin C=1.40 {best['margin_c_140']:+.6f}"
        )
    print()
    print("Largest MVDC block upper margins:")
    for row in sorted(block_rows, key=lambda item: item.q_upper_margin, reverse=True)[:8]:
        print(
            f"  pk={row.last_prime:>9} alpha={row.alpha:.2f} y={row.cutoff_y:>9} "
            f"Q={row.q_exact:+.6e} Q_upper_margin={row.q_upper_margin:+.3e} "
            f"max|u|={row.max_abs_u:.3e}"
        )


if __name__ == "__main__":
    main()
