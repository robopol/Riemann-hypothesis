from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_JSON = ROOT / "optimal_ledger_envelope_report.json"
DEFAULT_CSV = ROOT / "optimal_ledger_envelope_rows.csv"

EULER_GAMMA = 0.577215664901532860606512090082402431


@dataclass
class OptimalLedgerRow:
    x: int
    theta_x: float
    beta_error_e: float
    optimal_s: float
    optimal_log_n: float
    optimal_deficit_a: float
    optimal_bridge_b: float
    optimal_total_a_plus_b: float
    optimal_margin_vs_e: float
    active_prime_cutoff: float
    active_prime_count: int
    c_active_over_sqrt_x: float
    max_continuous_exponent: float
    floor_log_deficit_a: float
    floor_log_bridge_b: float
    floor_log_total_a_plus_b: float
    floor_log_margin_vs_e: float


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


def parse_int_list(value: str) -> list[int]:
    return [int(part.strip()) for part in value.split(",") if part.strip()]


def build_prefixes(primes: list[int]) -> tuple[dict[int, float], dict[int, float]]:
    log_beta = 0.0
    theta = 0.0
    log_beta_by_prime: dict[int, float] = {}
    theta_by_prime: dict[int, float] = {}
    for prime in primes:
        log_beta += -math.log1p(-1.0 / prime)
        theta += math.log(prime)
        log_beta_by_prime[prime] = log_beta
        theta_by_prime[prime] = theta
    return log_beta_by_prime, theta_by_prime


def build_prefix_arrays(primes: list[int]) -> tuple[list[float], list[float], list[float]]:
    theta = [0.0]
    f2 = [0.0]
    log_beta = [0.0]
    theta_total = 0.0
    f2_total = 0.0
    beta_total = 0.0
    for prime in primes:
        theta_total += math.log(prime)
        f2_total += -math.log1p(-1.0 / (prime * prime))
        beta_total += -math.log1p(-1.0 / prime)
        theta.append(theta_total)
        f2.append(f2_total)
        log_beta.append(beta_total)
    return theta, f2, log_beta


def beta_error(log_beta: float, x: int) -> float:
    return log_beta - EULER_GAMMA - math.log(math.log(x))


def bridge_reserve(log_n: float, x: int) -> float:
    if log_n <= x:
        return 0.0
    return math.log(math.log(log_n) / math.log(x))


def continuous_profile_values(
    primes: list[int],
    theta_prefix: list[float],
    f2_prefix: list[float],
    x_index: int,
    s: float,
) -> tuple[float, float, int, float]:
    """Evaluate a_p=max(1, s/log(p)-1) using prefix sums.

    For active primes p<=exp(s/2), p^{-(a_p+1)}=exp(-s), so the active
    deficit is constant per prime. The remaining primes have exponent 1.
    """
    active_cutoff = math.exp(0.5 * s) if s < 1400.0 else float("inf")
    active_count = min(x_index, bisect.bisect_right(primes, active_cutoff))

    theta_x = theta_prefix[x_index]
    theta_active = theta_prefix[active_count]
    log_n = theta_x + active_count * s - 2.0 * theta_active

    active_deficit = active_count * (-math.log1p(-math.exp(-s))) if active_count else 0.0
    tail_deficit = f2_prefix[x_index] - f2_prefix[active_count]
    deficit = active_deficit + tail_deficit
    max_exponent = s / math.log(2.0) - 1.0 if active_count else 1.0
    return log_n, deficit, active_count, max_exponent


def objective(
    primes: list[int],
    theta_prefix: list[float],
    f2_prefix: list[float],
    x_index: int,
    x: int,
    s: float,
) -> tuple[float, float, float, int, float]:
    log_n, deficit, active_count, max_exponent = continuous_profile_values(
        primes, theta_prefix, f2_prefix, x_index, s
    )
    if log_n < x:
        # The bridge theorem gives log(n)>x for a least counterexample.
        # Penalize profiles that do not reach the admissible region.
        return float("inf"), log_n, deficit, active_count, max_exponent
    bridge = bridge_reserve(log_n, x)
    return deficit + bridge, log_n, deficit, active_count, max_exponent


def find_optimal_s(
    primes: list[int],
    theta_prefix: list[float],
    f2_prefix: list[float],
    x_index: int,
    x: int,
) -> tuple[float, float, float, int, float]:
    lower = 2.0 * math.log(2.0)
    upper = 2.0 * math.log(x) + 8.0
    while objective(primes, theta_prefix, f2_prefix, x_index, x, upper)[1] < x:
        upper *= 1.5

    # Coarse scan brackets the best piecewise-smooth region.
    best_s = lower
    best_value = float("inf")
    scan_steps = 90
    for index in range(scan_steps + 1):
        s = lower + (upper - lower) * index / scan_steps
        value = objective(primes, theta_prefix, f2_prefix, x_index, x, s)[0]
        if value < best_value:
            best_value = value
            best_s = s

    step = (upper - lower) / scan_steps
    lo = max(lower, best_s - 2.0 * step)
    hi = min(upper, best_s + 2.0 * step)

    # Golden-section minimization.
    inv_phi = (math.sqrt(5.0) - 1.0) / 2.0
    inv_phi2 = (3.0 - math.sqrt(5.0)) / 2.0
    h = hi - lo
    c = lo + inv_phi2 * h
    d = lo + inv_phi * h
    fc = objective(primes, theta_prefix, f2_prefix, x_index, x, c)[0]
    fd = objective(primes, theta_prefix, f2_prefix, x_index, x, d)[0]
    for _ in range(80):
        if fc < fd:
            hi = d
            d = c
            fd = fc
            h = hi - lo
            c = lo + inv_phi2 * h
            fc = objective(primes, theta_prefix, f2_prefix, x_index, x, c)[0]
        else:
            lo = c
            c = d
            fc = fd
            h = hi - lo
            d = lo + inv_phi * h
            fd = objective(primes, theta_prefix, f2_prefix, x_index, x, d)[0]

    s = (lo + hi) / 2.0
    total, log_n, deficit, active_count, max_exponent = objective(
        primes, theta_prefix, f2_prefix, x_index, x, s
    )
    return s, log_n, deficit, active_count, max_exponent


def floor_log_profile(primes: list[int], x: int) -> tuple[float, float]:
    log_x = math.log(x)
    log_n = 0.0
    deficit = 0.0
    for prime in primes:
        if prime > x:
            break
        exponent = max(1, int(log_x // math.log(prime)))
        log_n += exponent * math.log(prime)
        deficit += -math.log1p(-(prime ** (-(exponent + 1))))
    return log_n, deficit


def build_rows(primes: list[int], targets: list[int]) -> list[OptimalLedgerRow]:
    theta_prefix, f2_prefix, log_beta_prefix = build_prefix_arrays(primes)
    rows: list[OptimalLedgerRow] = []
    for target in targets:
        x = nearest_prime_at_or_below(primes, target)
        x_index = bisect.bisect_right(primes, x)
        e_value = beta_error(log_beta_prefix[x_index], x)
        s, opt_log_n, opt_a, active_count, max_exponent = find_optimal_s(
            primes, theta_prefix, f2_prefix, x_index, x
        )
        opt_b = bridge_reserve(opt_log_n, x)
        opt_total = opt_a + opt_b

        fl_log_n, fl_a = floor_log_profile(primes, x)
        fl_b = bridge_reserve(fl_log_n, x)
        fl_total = fl_a + fl_b

        active_cutoff = math.exp(s / 2.0)
        rows.append(
            OptimalLedgerRow(
                x=x,
                theta_x=theta_prefix[x_index],
                beta_error_e=e_value,
                optimal_s=s,
                optimal_log_n=opt_log_n,
                optimal_deficit_a=opt_a,
                optimal_bridge_b=opt_b,
                optimal_total_a_plus_b=opt_total,
                optimal_margin_vs_e=opt_total - e_value,
                active_prime_cutoff=active_cutoff,
                active_prime_count=active_count,
                c_active_over_sqrt_x=active_cutoff / math.sqrt(x),
                max_continuous_exponent=max_exponent,
                floor_log_deficit_a=fl_a,
                floor_log_bridge_b=fl_b,
                floor_log_total_a_plus_b=fl_total,
                floor_log_margin_vs_e=fl_total - e_value,
            )
        )
    return rows


def write_csv(path: Path, rows: list[OptimalLedgerRow]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find the continuous lower envelope for A(n)+B_log(n) using all exponents."
    )
    parser.add_argument("--max-x", type=int, default=50_000_000)
    parser.add_argument("--targets", default="1000003,1999993,5000011,9999991,19999999,46909099,50000017")
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    targets = parse_int_list(args.targets)
    max_needed = max(args.max_x, max(targets))
    primes = primes_up_to(max_needed)
    rows = build_rows(primes, targets)

    payload = {
        "notes": [
            "This script optimizes the full ledger reserve A(n)+B_log(n), not a tail-only envelope.",
            "The continuous relaxed profile is a_p=max(1, s/log(p)-1), the water-filling solution for the deficit/log-size tradeoff.",
            "Because exponents are relaxed to real values, the resulting minimum is a lower bound for integer profiles in this model.",
            "This is still a model using full support p<=p_k and the bridge log(n)>=p_k; a proof must justify the candidate class assumptions.",
        ],
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    write_csv(args.csv, rows)

    print(f"primes up to {max_needed}: {len(primes)}")
    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print()
    print("Optimal continuous ledger envelope:")
    for row in rows:
        print(
            f"x={row.x:>9} E={row.beta_error_e:+.6e} "
            f"optA={row.optimal_deficit_a:.6e} optB={row.optimal_bridge_b:.6e} "
            f"optTotal={row.optimal_total_a_plus_b:.6e} "
            f"margin={row.optimal_margin_vs_e:+.6e} "
            f"cut/sqrt={row.c_active_over_sqrt_x:.4f} "
            f"floorMargin={row.floor_log_margin_vs_e:+.6e}"
        )


if __name__ == "__main__":
    main()
