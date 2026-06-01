from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_JSON = ROOT / "ca_cap_envelope_report.json"
DEFAULT_CSV = ROOT / "ca_cap_envelope_rows.csv"

EULER_GAMMA = 0.577215664901532860606512090082402431


@dataclass
class CACapEnvelopeRow:
    epsilon: float
    epsilon_low: float
    epsilon_high: float
    epsilon_in_support_interval: bool
    last_prime: int
    next_prime: int
    prime_count: int
    beta_error_e: float
    true_deficit_a: float
    true_bridge_b: float
    true_total_a_plus_b: float
    true_margin_vs_e: float
    cap_deficit_a: float
    cap_bridge_b: float
    cap_total_a_plus_b: float
    cap_margin_vs_e: float
    envelope_log_slack: float
    lower_log_n_slack: float
    max_upper_cap_minus_true: int
    max_true_minus_lower_cap: int
    upper_cap_violations: int
    lower_cap_violations: int
    first_upper_cap_excess_prime: int | None
    first_lower_cap_gap_prime: int | None
    max_true_exponent: int
    max_upper_cap_exponent: int
    max_lower_cap_exponent: int


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


def support_threshold(prime: int) -> float:
    """The epsilon threshold for a prime to have exponent at least one."""
    return math.log1p(1.0 / prime) / math.log(prime)


def ca_exponent(prime: int, epsilon: float) -> int:
    """Alaoglu-Erdos CA variational exponent for sigma(n)/n^(1+epsilon)."""
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


def bridge_reserve(log_n: float, x: int) -> float:
    if log_n <= x:
        return 0.0
    return math.log(math.log(log_n) / math.log(x))


def deficit_for_profile(profile: dict[int, int]) -> float:
    return sum(-math.log1p(-(prime ** (-(exponent + 1)))) for prime, exponent in profile.items())


def log_n_for_profile(profile: dict[int, int]) -> float:
    return sum(exponent * math.log(prime) for prime, exponent in profile.items())


def evaluate_caps(
    *,
    epsilon: float,
    profile: dict[int, int],
    support_primes: list[int],
    next_prime: int,
    log_beta_by_prime: dict[int, float],
) -> CACapEnvelopeRow:
    last_prime = support_primes[-1]
    epsilon_low = support_threshold(next_prime)
    epsilon_high = support_threshold(last_prime)
    epsilon_high_inside = math.nextafter(epsilon_high, 0.0)

    upper_caps = {prime: ca_exponent(prime, epsilon_low) for prime in support_primes}
    lower_caps = {
        prime: max(1, ca_exponent(prime, epsilon_high_inside))
        for prime in support_primes
    }

    true_log_n = log_n_for_profile(profile)
    true_a = deficit_for_profile(profile)
    true_b = bridge_reserve(true_log_n, last_prime)

    cap_a = deficit_for_profile(upper_caps)
    cap_log_n = log_n_for_profile(lower_caps)
    cap_b = bridge_reserve(cap_log_n, last_prime)

    error = beta_error(log_beta_by_prime[last_prime], last_prime)
    true_total = true_a + true_b
    cap_total = cap_a + cap_b

    upper_violations = 0
    lower_violations = 0
    max_upper_cap_minus_true = 0
    max_true_minus_lower_cap = 0
    first_upper_cap_excess_prime: int | None = None
    first_lower_cap_gap_prime: int | None = None
    for prime in support_primes:
        true_exp = profile[prime]
        upper_exp = upper_caps[prime]
        lower_exp = lower_caps[prime]
        if true_exp > upper_exp:
            upper_violations += 1
        if true_exp < lower_exp:
            lower_violations += 1
        upper_excess = upper_exp - true_exp
        lower_gap = true_exp - lower_exp
        if upper_excess > max_upper_cap_minus_true:
            max_upper_cap_minus_true = upper_excess
            first_upper_cap_excess_prime = prime
        if lower_gap > max_true_minus_lower_cap:
            max_true_minus_lower_cap = lower_gap
            first_lower_cap_gap_prime = prime

    return CACapEnvelopeRow(
        epsilon=epsilon,
        epsilon_low=epsilon_low,
        epsilon_high=epsilon_high,
        epsilon_in_support_interval=(epsilon_low <= epsilon < epsilon_high),
        last_prime=last_prime,
        next_prime=next_prime,
        prime_count=len(support_primes),
        beta_error_e=error,
        true_deficit_a=true_a,
        true_bridge_b=true_b,
        true_total_a_plus_b=true_total,
        true_margin_vs_e=true_total - max(error, 0.0),
        cap_deficit_a=cap_a,
        cap_bridge_b=cap_b,
        cap_total_a_plus_b=cap_total,
        cap_margin_vs_e=cap_total - max(error, 0.0),
        envelope_log_slack=true_a - cap_a,
        lower_log_n_slack=true_log_n - cap_log_n,
        max_upper_cap_minus_true=max_upper_cap_minus_true,
        max_true_minus_lower_cap=max_true_minus_lower_cap,
        upper_cap_violations=upper_violations,
        lower_cap_violations=lower_violations,
        first_upper_cap_excess_prime=first_upper_cap_excess_prime,
        first_lower_cap_gap_prime=first_lower_cap_gap_prime,
        max_true_exponent=max(profile.values()),
        max_upper_cap_exponent=max(upper_caps.values()),
        max_lower_cap_exponent=max(lower_caps.values()),
    )


def write_csv(path: Path, rows: list[CACapEnvelopeRow]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit certified CA exponent-cap envelopes for beta-envelope recovery."
    )
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
    rows: list[CACapEnvelopeRow] = []

    for epsilon in parse_float_list(args.epsilons):
        profile = ca_like_profile(primes, epsilon)
        if not profile:
            continue
        support_primes = sorted(profile)
        next_index = len(support_primes)
        if next_index >= len(primes):
            raise ValueError(
                "max-prime must include the full profile support and the next prime after it"
            )
        rows.append(
            evaluate_caps(
                epsilon=epsilon,
                profile=profile,
                support_primes=support_primes,
                next_prime=primes[next_index],
                log_beta_by_prime=log_beta_by_prime,
            )
        )

    min_cap_margin = min((row.cap_margin_vs_e for row in rows), default=float("nan"))
    max_identity_risk = max((row.upper_cap_violations + row.lower_cap_violations for row in rows), default=0)
    payload = {
        "notes": [
            "This is a numerical audit, not the analytic proof.",
            "Analytic proof target 1: if P(n)=x for a CA/SA candidate, epsilon lies in [tau(next_prime), tau(x)), tau(p)=log(1+1/p)/log(p).",
            "Analytic proof target 2: CA exponents a_p(epsilon) are non-increasing in epsilon.",
            "Then U_p=a_p(tau(next_prime)) is a certified upper exponent cap and gives A_cap <= A(n).",
            "The lower caps use the right endpoint from below plus contiguous support p<=x, giving L_cap <= log(n).",
            "A positive cap_margin_vs_e means this certified cap envelope numerically absorbs the measured sharp beta surplus E_+(x).",
        ],
        "min_cap_margin_vs_e": min_cap_margin,
        "max_cap_violation_count": max_identity_risk,
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    write_csv(args.csv, rows)

    print(f"primes up to {args.max_prime}: {len(primes)}")
    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print(f"min cap margin vs E_+: {min_cap_margin:+.6e}")
    print(f"max cap violation count: {max_identity_risk}")
    print()
    print("CA exponent-cap envelope:")
    for row in rows:
        print(
            f"eps={row.epsilon:g} pk={row.last_prime:>9} "
            f"E+={max(row.beta_error_e, 0.0):.6e} "
            f"A_cap={row.cap_deficit_a:.6e} B_cap={row.cap_bridge_b:.6e} "
            f"cap_margin={row.cap_margin_vs_e:+.6e} "
            f"slackA={row.envelope_log_slack:.3e} slackLogN={row.lower_log_n_slack:.3e} "
            f"ok_interval={row.epsilon_in_support_interval}"
        )


if __name__ == "__main__":
    main()
