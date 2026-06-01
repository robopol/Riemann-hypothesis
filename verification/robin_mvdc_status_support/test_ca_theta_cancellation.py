from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_JSON = ROOT / "ca_theta_cancellation_report.json"
DEFAULT_CSV = ROOT / "ca_theta_cancellation_rows.csv"

EULER_GAMMA = 0.577215664901532860606512090082402431


@dataclass
class CancellationRow:
    epsilon: float
    last_prime: int
    beta_error_e: float
    theta_minus_x: float
    ca_layer_h: float
    true_deficit_a: float
    bridge_b: float
    bridge_quadratic_remainder: float
    endpoint_term: float
    modified_error_w: float
    reserve_after_cancellation: float
    exact_gap: float
    reconstructed_gap: float
    reconstruction_error: float
    scaled_endpoint: float
    scaled_modified_error_w: float
    scaled_a_plus_h: float
    scaled_reserve_after_cancellation: float
    scaled_exact_gap: float


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


def ca_exponent(prime: int, epsilon: float) -> int:
    log_p = math.log(prime)
    numerator_log = (1.0 + epsilon) * log_p + math.log1p(-math.exp(-(1.0 + epsilon) * log_p))
    denominator_log = math.log(math.expm1(epsilon * log_p))
    exponent = math.floor((numerator_log - denominator_log) / log_p) - 1
    return max(0, exponent)


def ca_like_profile(all_primes: list[int], epsilon: float) -> dict[int, int]:
    profile: dict[int, int] = {}
    for prime in all_primes:
        exponent = ca_exponent(prime, epsilon)
        if exponent > 0:
            profile[prime] = exponent
        elif profile:
            break
    return profile


def build_prime_tables(primes: list[int]) -> tuple[dict[int, float], dict[int, float]]:
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


def evaluate_profile(
    epsilon: float,
    profile: dict[int, int],
    log_beta_by_prime: dict[int, float],
    theta_by_prime: dict[int, float],
) -> CancellationRow:
    primes = sorted(profile)
    x = primes[-1]
    log_x = math.log(x)
    log_n = 0.0
    deficit_a = 0.0
    for prime in primes:
        exponent = profile[prime]
        log_n += exponent * math.log(prime)
        deficit_a += -math.log1p(-(prime ** (-(exponent + 1))))

    theta_x = theta_by_prime[x]
    theta_minus_x = theta_x - x
    ca_layer_h = log_n - theta_x
    beta_error_e = log_beta_by_prime[x] - EULER_GAMMA - math.log(log_x)
    bridge_b = math.log(math.log(log_n) / log_x)

    endpoint_term = theta_minus_x / (x * log_x)
    modified_error_w = beta_error_e - endpoint_term

    linear_bridge = (theta_minus_x + ca_layer_h) / (x * log_x)
    bridge_quadratic_remainder = bridge_b - linear_bridge

    reserve_after_cancellation = deficit_a + ca_layer_h / (x * log_x) + bridge_quadratic_remainder
    exact_gap = deficit_a + bridge_b - beta_error_e
    reconstructed_gap = reserve_after_cancellation - modified_error_w
    scale = math.sqrt(x) * log_x

    return CancellationRow(
        epsilon=epsilon,
        last_prime=x,
        beta_error_e=beta_error_e,
        theta_minus_x=theta_minus_x,
        ca_layer_h=ca_layer_h,
        true_deficit_a=deficit_a,
        bridge_b=bridge_b,
        bridge_quadratic_remainder=bridge_quadratic_remainder,
        endpoint_term=endpoint_term,
        modified_error_w=modified_error_w,
        reserve_after_cancellation=reserve_after_cancellation,
        exact_gap=exact_gap,
        reconstructed_gap=reconstructed_gap,
        reconstruction_error=reconstructed_gap - exact_gap,
        scaled_endpoint=endpoint_term * scale,
        scaled_modified_error_w=modified_error_w * scale,
        scaled_a_plus_h=(deficit_a + ca_layer_h / (x * log_x)) * scale,
        scaled_reserve_after_cancellation=reserve_after_cancellation * scale,
        scaled_exact_gap=exact_gap * scale,
    )


def write_csv(path: Path, rows: list[CancellationRow]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Check the CA theta-endpoint cancellation identity.")
    parser.add_argument("--max-prime", type=int, default=50_000_000)
    parser.add_argument(
        "--epsilons",
        default="0.001,0.0005,0.0002,0.0001,0.00005,0.00002,0.00001,0.000005,0.000002,0.000001,0.0000005,0.0000002,0.0000001,0.00000005,0.00000002,0.00000001,0.000000005,0.000000002,0.000000001",
    )
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    primes = primes_up_to(args.max_prime)
    log_beta_by_prime, theta_by_prime = build_prime_tables(primes)

    rows: list[CancellationRow] = []
    for epsilon in parse_float_list(args.epsilons):
        profile = ca_like_profile(primes, epsilon)
        if profile:
            rows.append(evaluate_profile(epsilon, profile, log_beta_by_prime, theta_by_prime))

    payload = {
        "identity": "A+B_log-E = A+H/(x log x)+B2-(E-(theta(x)-x)/(x log x))",
        "notes": [
            "H is log(n_CA)-theta(x), the mass of the CA layers below the top layer.",
            "B2 is the exact nonlinear remainder after linearising the logarithmic bridge.",
            "reconstruction_error should be near floating-point roundoff.",
        ],
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    write_csv(args.csv, rows)

    print(f"primes up to {args.max_prime}: {len(primes)}")
    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print()
    print("CA theta cancellation:")
    for row in rows[-10:]:
        print(
            f"pk={row.last_prime:>9} "
            f"W={row.scaled_modified_error_w:6.3f} "
            f"A+H={row.scaled_a_plus_h:6.3f} "
            f"reserve={row.scaled_reserve_after_cancellation:6.3f} "
            f"gap={row.scaled_exact_gap:6.3f} "
            f"err={row.reconstruction_error:+.3e}"
        )


if __name__ == "__main__":
    main()
