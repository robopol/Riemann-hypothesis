from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
DEFAULT_SOURCE = ROOT / "source_copies" / "Riemann_test.py"
DEFAULT_OUTPUT = ROOT / "data" / "hcn_sa_bridge_test_results.json"

EULER_GAMMA = 0.5772156649015328606
E_GAMMA = math.exp(EULER_GAMMA)


def nth_prime_upper_bound(n: int) -> int:
    """Return a safe upper bound for the nth prime for the ranges used here."""
    if n < 6:
        return 15
    value = float(n)
    return int(math.ceil(value * (math.log(value) + math.log(math.log(value))) + 32))


def first_n_primes(n: int) -> list[int]:
    """Generate the first n primes using a bytearray sieve."""
    if n < 1:
        raise ValueError("n must be positive")

    limit = nth_prime_upper_bound(n)
    while True:
        sieve = bytearray(b"\x01") * (limit + 1)
        sieve[0:2] = b"\x00\x00"
        root = math.isqrt(limit)
        for value in range(2, root + 1):
            if sieve[value]:
                start = value * value
                sieve[start : limit + 1 : value] = b"\x00" * (((limit - start) // value) + 1)
        primes = [value for value, is_prime in enumerate(sieve) if is_prime]
        if len(primes) >= n:
            return primes[:n]
        limit *= 2


def checkpoint_counts(max_count: int, checkpoints: int) -> list[int]:
    """Build checkpoints with extra density near the values discussed in the audit notes."""
    linear = {max(1, round(max_count * index / checkpoints)) for index in range(1, checkpoints + 1)}
    fixed = {4_000, 40_000, 100_000, 300_000, 500_000}
    return sorted(count for count in linear | fixed if count <= max_count)


def local_sigma_over_n_factor(prime: int, exponent: int) -> float:
    """Return (1 + p + ... + p^a) / p^a without constructing p^a."""
    return prime / (prime - 1) - math.exp(-exponent * math.log(prime)) / (prime - 1)


def compute_log_profile(primes: list[int]) -> dict[str, float]:
    """Mirror Riemann_test.basic_sequence in log space.

    The original GUI constructs the huge candidate integer and its divisor sum.
    This version keeps the same greedy exponent walk but tracks log(N) and
    sigma(N)/N directly, so the million-prime bridge can be tested quickly.
    """
    log_n = 0.0
    log_sigma_over_n = 0.0
    log_beta = 0.0
    for prime in primes:
        log_prime = math.log(prime)
        log_n += log_prime
        log_sigma_over_n += math.log1p(1.0 / prime)
        log_beta += log_prime - math.log(prime - 1)

    sigma_over_n = math.exp(log_sigma_over_n)
    guy_robin = sigma_over_n / math.log(log_n)
    best_log_n = log_n
    best_sigma_over_n = sigma_over_n
    best_guy_robin = guy_robin

    consecutive_non_improvements = 0
    for prime in primes:
        exponent = 1
        while True:
            previous_factor = local_sigma_over_n_factor(prime, exponent)
            next_factor = local_sigma_over_n_factor(prime, exponent + 1)
            sigma_over_n *= next_factor / previous_factor
            log_n += math.log(prime)
            exponent += 1
            guy_robin = sigma_over_n / math.log(log_n)
            if guy_robin > best_guy_robin:
                consecutive_non_improvements = 0
                best_guy_robin = guy_robin
                best_log_n = log_n
                best_sigma_over_n = sigma_over_n
                continue
            consecutive_non_improvements += 1
            break
        if consecutive_non_improvements >= 2:
            break

    return {
        "log_n": best_log_n,
        "sigma_over_n": best_sigma_over_n,
        "guy_robin_index": best_guy_robin,
        "beta_k": math.exp(log_beta),
    }


def run_bridge_test(max_count: int, checkpoints: int, source: Path) -> dict[str, Any]:
    """Run the HCN/SA-like log-profile search and track both bridge conditions."""
    all_primes = first_n_primes(max_count)
    counts = checkpoint_counts(max_count, checkpoints)
    rows: list[dict[str, Any]] = []
    started = time.perf_counter()

    for step, current_count in enumerate(counts, start=1):
        step_started = time.perf_counter()
        primes = all_primes[:current_count]
        profile = compute_log_profile(primes)
        step_seconds = time.perf_counter() - step_started

        last_prime = primes[-1]
        log_last_prime = math.log(last_prime)
        log_n = float(profile["log_n"])
        sigma_over_n = float(profile["sigma_over_n"])
        guy_max = float(profile["guy_robin_index"])
        log_log_n = math.log(float(log_n))
        beta_k = float(profile["beta_k"])
        robin_rhs = E_GAMMA * log_log_n
        sigma_margin_to_robin = robin_rhs - float(sigma_over_n)
        beta_margin_to_robin = robin_rhs - beta_k
        bridge_margin = float(log_n) - last_prime
        delta_beta = beta_k / float(log_last_prime) - E_GAMMA
        epsilon_required_for_bridge = delta_beta / E_GAMMA
        required_bridge_power = 1.0 + epsilon_required_for_bridge
        strong_bridge_required_log_n = last_prime * math.exp(epsilon_required_for_bridge * float(log_last_prime))
        strong_bridge_room = float(log_n) - strong_bridge_required_log_n
        corrected_last_prime = strong_bridge_required_log_n
        corrected_last_prime_margin = strong_bridge_room

        rows.append(
            {
                "step": step,
                "prime_count": current_count,
                "last_prime": last_prime,
                "log_last_prime": float(log_last_prime),
                "log_n": float(log_n),
                "log_log_n": log_log_n,
                "n_digits": int(float(log_n) / math.log(10)) + 1,
                "bridge_log_n_minus_last_prime": bridge_margin,
                "bridge_holds_pk_lt_log_n": bridge_margin > 0.0,
                "last_prime_over_log_n": last_prime / float(log_n),
                "sigma_over_n": float(sigma_over_n),
                "guy_robin_index": float(guy_max),
                "e_gamma_margin": E_GAMMA - float(guy_max),
                "beta_k": beta_k,
                "beta_over_sigma": beta_k / float(sigma_over_n),
                "robin_rhs_sigma_scale": robin_rhs,
                "sigma_margin_to_robin_rhs": sigma_margin_to_robin,
                "beta_margin_to_robin_rhs": beta_margin_to_robin,
                "beta_sits_above_robin_rhs": beta_margin_to_robin < 0.0,
                "delta_beta": delta_beta,
                "epsilon_required_for_bridge": epsilon_required_for_bridge,
                "required_bridge_power": required_bridge_power,
                "strong_bridge_required_log_n": strong_bridge_required_log_n,
                "strong_bridge_room": strong_bridge_room,
                "strong_bridge_holds": strong_bridge_room >= 0.0,
                "corrected_last_prime": corrected_last_prime,
                "corrected_last_prime_margin": corrected_last_prime_margin,
                "corrected_last_prime_holds": corrected_last_prime_margin >= 0.0,
                "step_seconds": step_seconds,
            }
        )

        print(
            f"step={step:02d}/{len(counts)} count={current_count} "
            f"pk={last_prime} logN-pk={bridge_margin:+.3f} logN-Pcorr={corrected_last_prime_margin:+.3f} "
            f"G={float(guy_max):.10f} sigma_margin={sigma_margin_to_robin:+.6e} "
            f"digits={int(float(log_n) / math.log(10)) + 1} time={step_seconds:.2f}s",
            flush=True,
        )

    elapsed = time.perf_counter() - started
    max_guy_row = max(rows, key=lambda row: row["guy_robin_index"])
    min_bridge_row = min(rows, key=lambda row: row["bridge_log_n_minus_last_prime"])
    min_strong_bridge_row = min(rows, key=lambda row: row["strong_bridge_room"])
    return {
        "source_script": str(source.relative_to(ROOT)),
        "test_role": (
            "HCN/SA-like structured log-profile search. Unlike the last-prime beta test, this run "
            "tracks log(N), p_k, sigma(N)/N, the Guy Robin index, the bridge p_k < log(N), and "
            "the corrected last-prime bridge log(N) >= P_k^corr, where P_k^corr = p_k^(1 + delta_beta/e^gamma)."
        ),
        "implementation_note": (
            "This is a log-domain equivalent of the original Riemann_test.basic_sequence greedy exponent walk. "
            "It avoids constructing the huge integer N, but preserves the tracked log(N), sigma(N)/N, and G(N) values."
        ),
        "max_count": max_count,
        "checkpoints": len(rows),
        "requested_checkpoints": checkpoints,
        "euler_gamma": EULER_GAMMA,
        "e_gamma": E_GAMMA,
        "elapsed_seconds": elapsed,
        "bridge_failures": sum(1 for row in rows if not row["bridge_holds_pk_lt_log_n"]),
        "strong_bridge_failures": sum(1 for row in rows if not row["strong_bridge_holds"]),
        "corrected_last_prime_bridge_failures": sum(1 for row in rows if not row["corrected_last_prime_holds"]),
        "min_bridge_margin": min_bridge_row["bridge_log_n_minus_last_prime"],
        "min_bridge_margin_row": min_bridge_row,
        "min_strong_bridge_room": min_strong_bridge_row["strong_bridge_room"],
        "min_strong_bridge_room_row": min_strong_bridge_row,
        "min_corrected_last_prime_margin": min_strong_bridge_row["corrected_last_prime_margin"],
        "min_corrected_last_prime_margin_row": min_strong_bridge_row,
        "max_guy_robin_index": max_guy_row["guy_robin_index"],
        "min_margin_to_e_gamma": max_guy_row["e_gamma_margin"],
        "max_guy_row": max_guy_row,
        "exceedance_observed": any(row["guy_robin_index"] >= E_GAMMA for row in rows),
        "rows": rows,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the HCN/SA-like bridge tracking test.")
    parser.add_argument("--max-count", type=int, default=1_000_000)
    parser.add_argument("--checkpoints", type=int, default=24)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result = run_bridge_test(max_count=args.max_count, checkpoints=args.checkpoints, source=args.source)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(args.output)
    print(f"elapsed={result['elapsed_seconds']:.2f}s")
    print(f"bridge_failures={result['bridge_failures']}")
    print(f"max_guy_robin_index={result['max_guy_robin_index']:.10f}")


if __name__ == "__main__":
    main()
