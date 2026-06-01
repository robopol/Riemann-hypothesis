from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_JSON = ROOT / "low_prime_robin_region_report.json"
DEFAULT_CSV = ROOT / "low_prime_robin_region_worst.csv"

EULER_GAMMA = 0.577215664901532860606512090082402431
E_GAMMA = math.exp(EULER_GAMMA)
ROBIN_START = 5040


@dataclass(frozen=True)
class FrontierState:
    log_n: float
    log_sigma_over_n: float
    exponents: tuple[int, ...]


@dataclass
class WorstRow:
    rank: int
    n_decimal: str
    log_n: float
    largest_prime_factor: int
    sigma_over_n: float
    robin_boundary: float
    robin_index: float
    robin_margin: float
    exponents: str


def primes_below(limit: int) -> list[int]:
    primes: list[int] = []
    for candidate in range(2, limit):
        root = math.isqrt(candidate)
        if all(candidate % divisor for divisor in range(2, root + 1)):
            primes.append(candidate)
    return primes


def sigma_over_n_factor(prime: int, exponent: int) -> float:
    if exponent <= 0:
        return 1.0
    return (prime / (prime - 1.0)) * (1.0 - prime ** (-(exponent + 1)))


def beta_for_primes(primes: list[int]) -> float:
    beta = 1.0
    for prime in primes:
        beta *= prime / (prime - 1.0)
    return beta


def decimal_from_profile(primes: list[int], exponents: tuple[int, ...]) -> str:
    value = 1
    for prime, exponent in zip(primes, exponents):
        if exponent:
            value *= prime**exponent
    return str(value)


def largest_prime_factor(primes: list[int], exponents: tuple[int, ...]) -> int:
    for prime, exponent in reversed(list(zip(primes, exponents))):
        if exponent:
            return prime
    return 1


def format_profile(primes: list[int], exponents: tuple[int, ...]) -> str:
    parts = []
    for prime, exponent in zip(primes, exponents):
        if exponent:
            parts.append(f"{prime}^{exponent}")
    return " ".join(parts)


def prune_frontier(states: list[FrontierState], robin_log_start: float) -> list[FrontierState]:
    """Keep all below-threshold states and Pareto-prune only valid-domain states.

    A valid-domain state can safely dominate another valid-domain state if it has
    no larger log(n) and no smaller sigma(n)/n. States below 5040 are retained,
    because a below-domain state must not be allowed to erase a valid Robin
    candidate merely by having a smaller denominator.
    """

    small: list[FrontierState] = []
    large: list[FrontierState] = []
    for state in states:
        if state.log_n <= robin_log_start:
            small.append(state)
        else:
            large.append(state)

    large.sort(key=lambda state: (state.log_n, -state.log_sigma_over_n))
    pruned_large: list[FrontierState] = []
    best_log_sigma = -math.inf
    tolerance = 1e-15
    for state in large:
        if state.log_sigma_over_n > best_log_sigma + tolerance:
            pruned_large.append(state)
            best_log_sigma = state.log_sigma_over_n

    return small + pruned_large


def verify_low_prime_region(pk_threshold: int, top: int) -> dict:
    primes = primes_below(pk_threshold)
    beta_max = beta_for_primes(primes)

    # If log(n) is at least this value, then sigma(n)/n <= beta_max
    # is already below e^gamma log log(n). Only the finite region below
    # this logarithmic cutoff needs direct enumeration.
    log_n_cutoff = math.exp(beta_max / E_GAMMA)
    robin_log_start = math.log(ROBIN_START)

    states = [FrontierState(log_n=0.0, log_sigma_over_n=0.0, exponents=())]
    counts: list[dict] = []

    for prime in primes:
        log_prime = math.log(prime)
        next_states: list[FrontierState] = []
        for state in states:
            max_exponent = int((log_n_cutoff - state.log_n) // log_prime)
            for exponent in range(max_exponent + 1):
                new_log_n = state.log_n + exponent * log_prime
                if new_log_n > log_n_cutoff + 1e-12:
                    break
                new_log_sigma = state.log_sigma_over_n + math.log(
                    sigma_over_n_factor(prime, exponent)
                )
                next_states.append(
                    FrontierState(
                        log_n=new_log_n,
                        log_sigma_over_n=new_log_sigma,
                        exponents=state.exponents + (exponent,),
                    )
                )
        states = prune_frontier(next_states, robin_log_start)
        valid_count = sum(1 for state in states if state.log_n > robin_log_start)
        counts.append(
            {
                "prime": prime,
                "frontier_states": len(states),
                "valid_domain_states": valid_count,
            }
        )

    valid_states = [state for state in states if state.log_n > robin_log_start]
    rows: list[WorstRow] = []
    violations: list[WorstRow] = []

    scored: list[tuple[float, FrontierState]] = []
    for state in valid_states:
        sigma_over_n = math.exp(state.log_sigma_over_n)
        robin_boundary = E_GAMMA * math.log(state.log_n)
        robin_margin = robin_boundary - sigma_over_n
        scored.append((robin_margin, state))

    scored.sort(key=lambda item: item[0])

    for rank, (margin, state) in enumerate(scored[:top], start=1):
        sigma_over_n = math.exp(state.log_sigma_over_n)
        robin_boundary = E_GAMMA * math.log(state.log_n)
        row = WorstRow(
            rank=rank,
            n_decimal=decimal_from_profile(primes, state.exponents),
            log_n=state.log_n,
            largest_prime_factor=largest_prime_factor(primes, state.exponents),
            sigma_over_n=sigma_over_n,
            robin_boundary=robin_boundary,
            robin_index=sigma_over_n / math.log(state.log_n),
            robin_margin=margin,
            exponents=format_profile(primes, state.exponents),
        )
        rows.append(row)
        if margin <= 0.0:
            violations.append(row)

    return {
        "constants": {
            "euler_gamma": EULER_GAMMA,
            "e_gamma": E_GAMMA,
            "robin_start": ROBIN_START,
            "pk_threshold": pk_threshold,
            "largest_prime_checked": primes[-1] if primes else None,
            "prime_count": len(primes),
            "beta_max_for_low_prime_region": beta_max,
            "log_n_cutoff": log_n_cutoff,
        },
        "notes": [
            "The check covers all integers n > 5040 whose prime factors are all < pk_threshold.",
            "For log(n) >= log_n_cutoff, the universal bound sigma(n)/n <= beta_max already proves Robin.",
            "For the remaining finite region, a Pareto frontier over exact prime-exponent profiles is enumerated.",
            "The frontier keeps all below-domain states and prunes only valid-domain states dominated in both log(n) and sigma(n)/n.",
        ],
        "frontier_counts": counts,
        "summary": {
            "valid_frontier_states": len(valid_states),
            "violations": len([item for item in scored if item[0] <= 0.0]),
            "worst_margin": scored[0][0] if scored else None,
            "worst_n": decimal_from_profile(primes, scored[0][1].exponents) if scored else None,
            "worst_largest_prime_factor": largest_prime_factor(primes, scored[0][1].exponents)
            if scored
            else None,
        },
        "worst_rows": [asdict(row) for row in rows],
        "violations_sample": [asdict(row) for row in violations[:top]],
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify Robin's inequality in the low-prime region p_k < 286."
    )
    parser.add_argument("--pk-threshold", type=int, default=286)
    parser.add_argument("--top", type=int, default=25)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = verify_low_prime_region(args.pk_threshold, args.top)
    args.json.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    write_csv(args.csv, result["worst_rows"])

    constants = result["constants"]
    summary = result["summary"]
    print(f"checked p < {constants['pk_threshold']} ({constants['prime_count']} primes)")
    print(f"beta_max={constants['beta_max_for_low_prime_region']:.15f}")
    print(f"log_n_cutoff={constants['log_n_cutoff']:.15f}")
    print(f"valid_frontier_states={summary['valid_frontier_states']}")
    print(f"violations={summary['violations']}")
    print(f"worst_margin={summary['worst_margin']:.15e}")
    print(f"worst_n={summary['worst_n']}")
    print(f"worst_largest_prime_factor={summary['worst_largest_prime_factor']}")
    print(args.json)
    print(args.csv)


if __name__ == "__main__":
    main()
