from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_JSON = ROOT / "ca_like_ledger_report.json"
DEFAULT_CSV = ROOT / "ca_like_ledger_rows.csv"

EULER_GAMMA = 0.577215664901532860606512090082402431


@dataclass
class CALikeLedgerRow:
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
    first_unit_prime: int
    first_unit_over_sqrt_last: float
    last_nonunit_prime: int | None
    max_exponent: int
    log_n_minus_last_prime: float


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


def ca_like_profile(all_primes: list[int], epsilon: float) -> dict[int, int]:
    """Build the CA variational profile using the Alaoglu-Erdos exponent formula."""
    profile: dict[int, int] = {}
    for prime in all_primes:
        exponent = ca_exponent(prime, epsilon)
        if exponent > 0:
            profile[prime] = exponent
        elif profile:
            break
    return profile


def evaluate_profile(epsilon: float, profile: dict[int, int], log_beta_by_prime: dict[int, float]) -> CALikeLedgerRow:
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
        raise ValueError("Profile has no unit tail.")
    bridge = math.log(math.log(log_n) / math.log(last_prime)) if log_n > last_prime else 0.0
    error = beta_error(log_beta_by_prime[last_prime], last_prime)
    scale = math.sqrt(last_prime) * math.log(last_prime)
    total = deficit + bridge
    return CALikeLedgerRow(
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
        first_unit_prime=first_unit,
        first_unit_over_sqrt_last=first_unit / math.sqrt(last_prime),
        last_nonunit_prime=last_nonunit,
        max_exponent=max_exponent,
        log_n_minus_last_prime=log_n - last_prime,
    )


def write_csv(path: Path, rows: list[CALikeLedgerRow]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate the exact ledger on CA-like variational profiles.")
    parser.add_argument("--max-prime", type=int, default=80_000_000)
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
    log_beta_by_prime = build_log_beta_by_prime(primes)
    rows: list[CALikeLedgerRow] = []
    for epsilon in parse_float_list(args.epsilons):
        profile = ca_like_profile(primes, epsilon)
        if profile:
            rows.append(evaluate_profile(epsilon, profile, log_beta_by_prime))

    payload = {
        "notes": [
            "Profiles are local CA-like maximizers of sigma(n)/n^(1+epsilon).",
            "This evaluates the exact ledger E(p_k) <= A(n)+B_log(n) on those profiles.",
            "A positive margin_vs_e means the measured beta error is absorbed by the actual profile reserve.",
            "This is numerical evidence for the structural envelope, not the analytic proof.",
        ],
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    write_csv(args.csv, rows)

    print(f"primes up to {args.max_prime}: {len(primes)}")
    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print()
    print("CA-like exact ledger:")
    for row in rows:
        print(
            f"eps={row.epsilon:g} pk={row.last_prime:>9} "
            f"E={row.beta_error_e:+.6e} A={row.true_deficit_a:.6e} B={row.bridge_b:.6e} "
            f"A+B={row.total_a_plus_b:.6e} margin={row.margin_vs_e:+.6e} "
            f"C_E={row.c_beta_actual:.4f} C_res={row.c_reserve_actual:.4f} "
            f"j1/sqrt={row.first_unit_over_sqrt_last:.4f}"
        )


if __name__ == "__main__":
    main()
