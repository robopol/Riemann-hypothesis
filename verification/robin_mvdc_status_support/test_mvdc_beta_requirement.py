from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_JSON = ROOT / "mvdc_beta_requirement_report.json"
DEFAULT_CSV = ROOT / "mvdc_beta_requirement_profiles.csv"
DEFAULT_BLOCK_CSV = ROOT / "mvdc_beta_requirement_blocks.csv"

EULER_GAMMA = 0.577215664901532860606512090082402431
E_GAMMA = math.exp(EULER_GAMMA)


@dataclass
class BetaPoint:
    last_prime: int
    prime_count: int
    log_pk: float
    log_beta: float
    beta_error_E: float
    beta_required_eta: float
    rosser_eta_half: float
    rosser_eta_one: float
    required_over_rosser_half: float
    beta_over_egamma_log: float


@dataclass
class ProfileRow:
    family: str
    label: str
    last_prime: int
    prime_count: int
    log_n: float
    bridge_log_n_minus_pk: float
    max_exponent: int
    j1_count: int
    beta_error_E: float
    beta_required_eta: float
    true_deficit_A: float
    rosser_eta_half: float
    rosser_eta_one: float
    mvdc_requirement_margin: float
    rosser_half_margin: float
    rosser_one_margin: float
    sigma_over_n_over_egamma_log_pk: float
    robin_margin: float


@dataclass
class BlockRow:
    last_prime: int
    cutoff_y: int
    alpha: float
    block_prime_count: int
    beta_error_E_x: float
    beta_error_E_y: float
    q_exact: float
    decomposition_error: float
    raw_block_log_product: float
    main_loglog_increment: float
    block_mean_log_factor: float
    block_second_central_moment: float


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


def parse_float_list(value: str) -> list[float]:
    return [float(part.strip()) for part in value.split(",") if part.strip()]


def build_log_beta_by_prime(primes: list[int]) -> dict[int, float]:
    log_beta = 0.0
    values: dict[int, float] = {}
    for prime in primes:
        log_beta += -math.log1p(-1.0 / prime)
        values[prime] = log_beta
    return values


def beta_error(log_beta: float, x: int) -> float:
    return log_beta - EULER_GAMMA - math.log(math.log(x))


def rosser_eta(x: int, c_m: float) -> float:
    log_x = math.log(x)
    return math.log1p(c_m / (log_x * log_x))


def profile_floor_log(primes: list[int], last_prime: int, scale: float = 1.0) -> dict[int, int]:
    log_pk = math.log(last_prime)
    return {
        prime: max(1, int(math.floor(scale * log_pk / math.log(prime))))
        for prime in primes
    }


def profile_all_j(primes: list[int], exponent: int) -> dict[int, int]:
    return {prime: exponent for prime in primes}


def profile_primorial(primes: list[int]) -> dict[int, int]:
    return {prime: 1 for prime in primes}


def ca_like_profile(all_primes: list[int], epsilon: float) -> dict[int, int]:
    """Build a local CA-like variational profile for sigma(n)/n^(1+epsilon)."""
    profile: dict[int, int] = {}
    for prime in all_primes:
        exponent = 0
        while True:
            numerator = prime ** (exponent + 2) - 1
            denominator = prime ** (exponent + 1) - 1
            gain = (numerator / denominator) / (prime ** (1.0 + epsilon))
            if gain <= 1.0:
                break
            exponent += 1
        if exponent > 0:
            profile[prime] = exponent
        elif profile:
            break
    return profile


def evaluate_beta_points(primes: list[int], log_beta_by_prime: dict[int, float], targets: list[int]) -> list[BetaPoint]:
    rows: list[BetaPoint] = []
    for target in targets:
        last_prime = nearest_prime_at_or_below(primes, target)
        index = primes.index(last_prime)
        log_pk = math.log(last_prime)
        log_beta = log_beta_by_prime[last_prime]
        error = beta_error(log_beta, last_prime)
        required = max(0.0, error)
        eta_half = rosser_eta(last_prime, 0.5)
        eta_one = rosser_eta(last_prime, 1.0)
        rows.append(
            BetaPoint(
                last_prime=last_prime,
                prime_count=index + 1,
                log_pk=log_pk,
                log_beta=log_beta,
                beta_error_E=error,
                beta_required_eta=required,
                rosser_eta_half=eta_half,
                rosser_eta_one=eta_one,
                required_over_rosser_half=required / eta_half if eta_half > 0 else float("nan"),
                beta_over_egamma_log=math.exp(error),
            )
        )
    return rows


def evaluate_profile(
    family: str,
    label: str,
    profile: dict[int, int],
    log_beta_by_prime: dict[int, float],
) -> ProfileRow:
    primes = sorted(profile)
    last_prime = primes[-1]
    log_pk = math.log(last_prime)
    log_beta = log_beta_by_prime[last_prime]
    error = beta_error(log_beta, last_prime)
    required = max(0.0, error)
    log_deficit = 0.0
    log_n = 0.0
    max_exponent = 0
    j1_count = 0

    for prime in primes:
        exponent = profile[prime]
        max_exponent = max(max_exponent, exponent)
        j1_count += 1 if exponent == 1 else 0
        log_n += exponent * math.log(prime)
        log_deficit += math.log1p(-(prime ** (-(exponent + 1))))

    true_deficit = -log_deficit
    eta_half = rosser_eta(last_prime, 0.5)
    eta_one = rosser_eta(last_prime, 1.0)
    sigma_over_egamma_log_pk = math.exp(error - true_deficit)
    robin_margin = (
        E_GAMMA * math.log(log_n) - math.exp(log_beta + log_deficit)
        if log_n > 1.0
        else float("nan")
    )

    return ProfileRow(
        family=family,
        label=label,
        last_prime=last_prime,
        prime_count=len(primes),
        log_n=log_n,
        bridge_log_n_minus_pk=log_n - last_prime,
        max_exponent=max_exponent,
        j1_count=j1_count,
        beta_error_E=error,
        beta_required_eta=required,
        true_deficit_A=true_deficit,
        rosser_eta_half=eta_half,
        rosser_eta_one=eta_one,
        mvdc_requirement_margin=true_deficit - required,
        rosser_half_margin=true_deficit - eta_half,
        rosser_one_margin=true_deficit - eta_one,
        sigma_over_n_over_egamma_log_pk=sigma_over_egamma_log_pk,
        robin_margin=robin_margin,
    )


def build_profile_rows(
    primes: list[int],
    log_beta_by_prime: dict[int, float],
    targets: list[int],
    ca_epsilons: list[float],
    scaled_profiles: list[float],
) -> list[ProfileRow]:
    rows: list[ProfileRow] = []
    for target in targets:
        last_prime = nearest_prime_at_or_below(primes, target)
        support = [prime for prime in primes if prime <= last_prime]
        rows.append(evaluate_profile("floor_log_lcm", f"pk={last_prime}", profile_floor_log(support, last_prime), log_beta_by_prime))
        rows.append(evaluate_profile("primorial", f"pk={last_prime}", profile_primorial(support), log_beta_by_prime))
        rows.append(evaluate_profile("all_j2", f"pk={last_prime}", profile_all_j(support, 2), log_beta_by_prime))
        for scale in scaled_profiles:
            rows.append(
                evaluate_profile(
                    f"scaled_floor_log_{scale:g}",
                    f"pk={last_prime}",
                    profile_floor_log(support, last_prime, scale),
                    log_beta_by_prime,
                )
            )

    for epsilon in ca_epsilons:
        profile = ca_like_profile(primes, epsilon)
        if profile:
            rows.append(evaluate_profile("ca_like", f"eps={epsilon:g}", profile, log_beta_by_prime))
    return rows


def evaluate_blocks(
    primes: list[int],
    log_beta_by_prime: dict[int, float],
    targets: list[int],
    alphas: list[float],
) -> list[BlockRow]:
    rows: list[BlockRow] = []
    for target in targets:
        x = nearest_prime_at_or_below(primes, target)
        e_x = beta_error(log_beta_by_prime[x], x)
        for alpha in alphas:
            y_target = int(x ** alpha)
            y = nearest_prime_at_or_below(primes, y_target)
            if y >= x:
                continue
            block = [prime for prime in primes if y < prime <= x]
            if not block:
                continue
            e_y = beta_error(log_beta_by_prime[y], y)
            raw = sum(-math.log1p(-1.0 / prime) for prime in block)
            main = math.log(math.log(x) / math.log(y))
            q_exact = raw - main
            values = [-math.log1p(-1.0 / prime) for prime in block]
            mean = sum(values) / len(values)
            second = sum((value - mean) ** 2 for value in values)
            rows.append(
                BlockRow(
                    last_prime=x,
                    cutoff_y=y,
                    alpha=alpha,
                    block_prime_count=len(block),
                    beta_error_E_x=e_x,
                    beta_error_E_y=e_y,
                    q_exact=q_exact,
                    decomposition_error=(e_y + q_exact) - e_x,
                    raw_block_log_product=raw,
                    main_loglog_increment=main,
                    block_mean_log_factor=mean,
                    block_second_central_moment=second,
                )
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
        description="Numerically audit the sharp beta-error that an MVDC replacement for Rosser would need."
    )
    parser.add_argument("--max-pk", type=int, default=2_000_000)
    parser.add_argument("--targets", default="283,1009,10007,100003,199999,1000003,1999993")
    parser.add_argument("--ca-epsilons", default="0.001,0.0005,0.0002,0.0001,0.00005,0.00002,0.00001,0.000005,0.000002,0.000001")
    parser.add_argument("--scaled-profiles", default="")
    parser.add_argument("--block-alphas", default="0.5,0.75,0.9")
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--block-csv", type=Path, default=DEFAULT_BLOCK_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    targets = parse_int_list(args.targets)
    max_needed = max([args.max_pk, *targets])
    primes = primes_up_to(max_needed)
    log_beta_by_prime = build_log_beta_by_prime(primes)

    beta_points = evaluate_beta_points(primes, log_beta_by_prime, targets)
    profile_rows = build_profile_rows(
        primes=primes,
        log_beta_by_prime=log_beta_by_prime,
        targets=targets,
        ca_epsilons=parse_float_list(args.ca_epsilons),
        scaled_profiles=parse_float_list(args.scaled_profiles),
    )
    block_rows = evaluate_blocks(
        primes=primes,
        log_beta_by_prime=log_beta_by_prime,
        targets=targets,
        alphas=parse_float_list(args.block_alphas),
    )

    payload = {
        "constants": {
            "euler_gamma": EULER_GAMMA,
            "e_gamma": E_GAMMA,
        },
        "notes": [
            "E(x)=log(beta(x))-gamma-log(log(x)).",
            "beta_required_eta=max(E(x),0). This is the sharp beta-error MVDC must upper-bound for compensation.",
            "true_deficit_A=-log(D(n)) is the maximum envelope deficit available without dropping below sigma(n)/n.",
            "mvdc_requirement_margin=true_deficit_A-beta_required_eta; nonnegative means a sharp beta bound would numerically suffice.",
            "rosser_half_margin=true_deficit_A-log(1+(1/2)/(log x)^2); negative confirms the Rosser surplus is too coarse for that profile.",
            "Q(Y,x) is the exact renormalised finite block: sum_{Y<p<=x}-log(1-1/p)-log(log x/log Y).",
        ],
        "args": {
            "max_pk": args.max_pk,
            "targets": targets,
            "ca_epsilons": parse_float_list(args.ca_epsilons),
            "scaled_profiles": parse_float_list(args.scaled_profiles),
            "block_alphas": parse_float_list(args.block_alphas),
        },
        "beta_points": [asdict(row) for row in beta_points],
        "profile_rows": [asdict(row) for row in profile_rows],
        "block_rows": [asdict(row) for row in block_rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    write_csv(args.csv, profile_rows)
    write_csv(args.block_csv, block_rows)

    print(f"primes up to {max_needed}: {len(primes)}")
    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print(f"wrote {args.block_csv}")
    print()
    print("selected beta errors:")
    for row in beta_points:
        print(
            f"pk={row.last_prime:>8} E={row.beta_error_E:+.6e} "
            f"E+={row.beta_required_eta:.6e} eta_R(1/2)={row.rosser_eta_half:.6e} "
            f"ratio={row.required_over_rosser_half:.4f}"
        )
    print()
    print("critical profile margins:")
    interesting = [
        row for row in profile_rows
        if row.family in {"floor_log_lcm", "ca_like"}
    ]
    for row in interesting:
        print(
            f"{row.family:13s} {row.label:12s} pk={row.last_prime:>8} "
            f"A={row.true_deficit_A:.6e} E+={row.beta_required_eta:.6e} "
            f"MVDC_margin={row.mvdc_requirement_margin:+.6e} "
            f"Rosser_half_margin={row.rosser_half_margin:+.6e}"
        )


if __name__ == "__main__":
    main()
