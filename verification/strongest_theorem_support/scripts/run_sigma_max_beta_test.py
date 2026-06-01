from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
DEFAULT_SOURCE = ROOT / "source_copies" / "sigma-max_test_robopol_raw.py"
DEFAULT_OUTPUT = ROOT / "data" / "sigma_max_beta_results.json"

EULER_GAMMA = 0.5772156649015328606
E_GAMMA = math.exp(EULER_GAMMA)

# The downloaded GUI script uses this value under the name e_gama. It is the
# Robin index at 5040, not Euler's e^gamma, so the adapted test keeps it separate.
LEGACY_ROBIN_5040_INDEX = 1.7909733665348811334
DELTA_TO_TRUE_E_GAMMA = LEGACY_ROBIN_5040_INDEX - E_GAMMA


def nth_prime_upper_bound(n: int) -> int:
    """Return a safe upper bound for the nth prime for the ranges used here."""
    if n < 6:
        return 15
    n_float = float(n)
    return int(math.ceil(n_float * (math.log(n_float) + math.log(math.log(n_float))) + 32))


def first_n_primes(n: int) -> list[int]:
    """Generate the first n primes with a bytearray sieve."""
    if n < 1:
        raise ValueError("n must be positive")

    limit = nth_prime_upper_bound(n)
    while True:
        sieve = bytearray(b"\x01") * (limit + 1)
        sieve[0:2] = b"\x00\x00"
        root = int(math.isqrt(limit))
        for value in range(2, root + 1):
            if sieve[value]:
                start = value * value
                sieve[start : limit + 1 : value] = b"\x00" * (((limit - start) // value) + 1)
        primes = [value for value, is_prime in enumerate(sieve) if is_prime]
        if len(primes) >= n:
            return primes[:n]
        limit *= 2


def checkpoint_counts(max_count: int, checkpoints: int) -> list[int]:
    """Build useful checkpoints, including the theorem thresholds mentioned in the notes."""
    linear = {max(1, round(max_count * index / checkpoints)) for index in range(1, checkpoints + 1)}
    fixed = {10, 100, 1_000, 10_000, 100_000, 300_000, 500_000}
    return sorted(count for count in linear | fixed if count <= max_count)


def run_beta_test(max_count: int, checkpoints: int) -> dict[str, Any]:
    """Compute beta_k and the last-prime bounds without constructing huge integers."""
    started = time.perf_counter()
    primes = first_n_primes(max_count)
    targets = checkpoint_counts(max_count, checkpoints)
    target_set = set(targets)

    rows: list[dict[str, Any]] = []
    log_beta = 0.0
    theta = 0.0

    for index, prime in enumerate(primes, start=1):
        log_beta += math.log(prime) - math.log(prime - 1)
        theta += math.log(prime)
        if index not in target_set:
            continue

        beta_value = math.exp(log_beta)
        last_prime_log = math.log(prime)
        e_gamma_bound = E_GAMMA * last_prime_log
        legacy_bound = LEGACY_ROBIN_5040_INDEX * last_prime_log
        delta_required_for_beta = beta_value / last_prime_log - E_GAMMA
        epsilon_required_for_bridge = delta_required_for_beta / E_GAMMA
        rows.append(
            {
                "prime_count": index,
                "last_prime": prime,
                "log_last_prime": last_prime_log,
                "theta_last_prime": theta,
                "primorial_digit_estimate": int(theta / math.log(10)) + 1,
                "beta": beta_value,
                "e_gamma_log_last_prime": e_gamma_bound,
                "delta_corrected_log_last_prime": legacy_bound,
                "legacy_constant_log_last_prime": legacy_bound,
                "margin_to_true_e_gamma_bound": e_gamma_bound - beta_value,
                "margin_to_delta_corrected_bound": legacy_bound - beta_value,
                "margin_to_legacy_bound": legacy_bound - beta_value,
                "delta_required_for_beta": delta_required_for_beta,
                "epsilon_required_for_bridge": epsilon_required_for_bridge,
                "delta0_room_over_required_delta": DELTA_TO_TRUE_E_GAMMA - delta_required_for_beta,
                "delta0_bridge_power": 1.0 + DELTA_TO_TRUE_E_GAMMA / E_GAMMA,
                "required_bridge_power": 1.0 + epsilon_required_for_bridge,
                "beta_over_true_e_gamma_bound": beta_value / e_gamma_bound,
                "beta_over_delta_corrected_bound": beta_value / legacy_bound,
                "beta_over_legacy_bound": beta_value / legacy_bound,
            }
        )

    elapsed = time.perf_counter() - started
    theorem_rows = [row for row in rows if row["prime_count"] >= 100]
    closest_true = min(rows, key=lambda row: row["margin_to_true_e_gamma_bound"])
    closest_legacy = min(rows, key=lambda row: row["margin_to_legacy_bound"])
    closest_true_theorem_range = min(theorem_rows, key=lambda row: row["margin_to_true_e_gamma_bound"])
    closest_legacy_theorem_range = min(theorem_rows, key=lambda row: row["margin_to_legacy_bound"])
    return {
        "source_script": str(DEFAULT_SOURCE.relative_to(ROOT)),
        "adaptation_note": (
            "Headless numerical adaptation of the downloaded GUI sigma-max_test.py. "
            "It computes beta_k by the last prime only and keeps Euler's true e^gamma "
            "separate from the larger legacy GUI constant."
        ),
        "max_count": max_count,
        "checkpoints": len(rows),
        "requested_checkpoints": checkpoints,
        "elapsed_seconds": elapsed,
        "euler_gamma": EULER_GAMMA,
        "e_gamma": E_GAMMA,
        "delta_corrected_constant": LEGACY_ROBIN_5040_INDEX,
        "delta_to_true_e_gamma": DELTA_TO_TRUE_E_GAMMA,
        "deficit_needed_to_remove_delta": E_GAMMA / LEGACY_ROBIN_5040_INDEX,
        "true_e_gamma_exceedance_observed": any(row["margin_to_true_e_gamma_bound"] < 0 for row in rows),
        "legacy_exceedance_observed": any(row["margin_to_legacy_bound"] < 0 for row in rows),
        "theorem_range_start_prime_count": 100,
        "true_e_gamma_exceedance_observed_from_100": any(
            row["margin_to_true_e_gamma_bound"] < 0 for row in theorem_rows
        ),
        "legacy_exceedance_observed_from_100": any(row["margin_to_legacy_bound"] < 0 for row in theorem_rows),
        "closest_true_e_gamma_row": closest_true,
        "closest_legacy_row": closest_legacy,
        "closest_true_e_gamma_row_from_100": closest_true_theorem_range,
        "closest_legacy_row_from_100": closest_legacy_theorem_range,
        "proof_dependency_check": {
            "consolidated_proof_uses_appendix_1_bridge": True,
            "bridge": "p_k < log n for the structured candidates",
            "appendix_2_role": "Mertens bound plus the multiplicative deficit/tail mechanism",
        },
        "rows": rows,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the beta-envelope test headlessly.")
    parser.add_argument("--max-count", type=int, default=500_000)
    parser.add_argument("--checkpoints", type=int, default=12)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = run_beta_test(max_count=args.max_count, checkpoints=args.checkpoints)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(args.output)
    print(f"elapsed={result['elapsed_seconds']:.2f}s")
    print(f"true_e_gamma_exceedance_observed={result['true_e_gamma_exceedance_observed']}")
    print(f"legacy_exceedance_observed={result['legacy_exceedance_observed']}")
    print(f"true_e_gamma_exceedance_observed_from_100={result['true_e_gamma_exceedance_observed_from_100']}")
    print(f"legacy_exceedance_observed_from_100={result['legacy_exceedance_observed_from_100']}")


if __name__ == "__main__":
    main()
