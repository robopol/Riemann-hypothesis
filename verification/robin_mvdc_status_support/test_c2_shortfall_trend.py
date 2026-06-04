from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
from array import array
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_JSON = ROOT / "c2_shortfall_trend_report.json"
DEFAULT_CSV = ROOT / "c2_shortfall_trend_rows.csv"

EULER_GAMMA = 0.577215664901532860606512090082402431
MEISSEL_MERTENS_PRIME_CONSTANT = 0.261497212847642783755426838608695859


@dataclass
class CATrendPoint:
    epsilon: float
    last_prime: int
    next_prime: int
    prime_count: int
    beta_error_e: float
    reserve_a_plus_b: float
    log_n_minus_last_prime: float
    max_exponent: int


@dataclass
class C2ShortfallTrendRow:
    y_prime: int
    x_prime: int
    block_prime_count: int
    previous_upper_e_y: float
    reserve_x: float
    moment_1: float
    m1_step_threshold: float
    required_constant: float
    actual_constant: float
    finite_c2_constant: float
    c2_shortfall: float
    additive_shortfall: float
    dusart_2018_constant: float
    axler_2018_constant: float | None
    axler_2018_margin: float | None
    actual_margin: float
    cumulative_upper_e_x: float
    beta_error_e_x: float
    upper_minus_actual_e: float


def parse_float_list(value: str) -> list[float]:
    return [float(part.strip()) for part in value.split(",") if part.strip()]


def primes_up_to(limit: int) -> array:
    if limit < 2:
        return array("I")
    size = limit // 2 + 1
    sieve = bytearray(b"\x01") * size
    sieve[0] = 0
    root = math.isqrt(limit)
    for odd in range(3, root + 1, 2):
        index = odd // 2
        if sieve[index]:
            start = (odd * odd) // 2
            step = odd
            sieve[start::step] = b"\x00" * (((size - 1 - start) // step) + 1)
    primes = array("I", [2])
    primes.extend((2 * index + 1 for index in range(1, size) if sieve[index]))
    return primes


def ca_exponent(prime: int, epsilon: float) -> int:
    log_p = math.log(prime)
    numerator_log = (1.0 + epsilon) * log_p + math.log1p(
        -math.exp(-(1.0 + epsilon) * log_p)
    )
    denominator_log = math.log(math.expm1(epsilon * log_p))
    exponent = math.floor((numerator_log - denominator_log) / log_p) - 1
    return max(0, exponent)


def beta_error(log_beta: float, x: int) -> float:
    return log_beta - EULER_GAMMA - math.log(math.log(x))


def build_prefixes(primes: array) -> tuple[array, array, array, array]:
    prefix_log_beta = array("d", [0.0])
    prefix_inv_p = array("d", [0.0])
    prefix_c2 = array("d", [0.0])
    prefix_s_minus1 = array("d", [0.0])
    log_beta = 0.0
    inv_sum = 0.0
    c2_sum = 0.0
    s_minus1_sum = 0.0
    for prime in primes:
        log_beta += -math.log1p(-1.0 / prime)
        inv_sum += 1.0 / prime
        c2_sum += 1.0 / (prime * (prime - 1.0))
        s_minus1_sum += 1.0 / (prime - 1.0)
        prefix_log_beta.append(log_beta)
        prefix_inv_p.append(inv_sum)
        prefix_c2.append(c2_sum)
        prefix_s_minus1.append(s_minus1_sum)
    return prefix_log_beta, prefix_inv_p, prefix_c2, prefix_s_minus1


def evaluate_ca_point(primes: array, prefix_log_beta: array, epsilon: float) -> CATrendPoint | None:
    log_n = 0.0
    deficit = 0.0
    max_exponent = 0
    last_prime = 0
    next_prime = 0
    prime_count = 0
    seen_support = False
    for prime in primes:
        exponent = ca_exponent(prime, epsilon)
        if exponent > 0:
            seen_support = True
            last_prime = prime
            prime_count += 1
            max_exponent = max(max_exponent, exponent)
            log_n += exponent * math.log(prime)
            deficit += -math.log1p(-(prime ** (-(exponent + 1))))
        elif seen_support:
            next_prime = prime
            break
    if not seen_support or next_prime == 0:
        return None

    index = bisect.bisect_right(primes, last_prime)
    error = beta_error(prefix_log_beta[index], last_prime)
    bridge = math.log(math.log(log_n) / math.log(last_prime)) if log_n > last_prime else 0.0
    return CATrendPoint(
        epsilon=epsilon,
        last_prime=last_prime,
        next_prime=next_prime,
        prime_count=prime_count,
        beta_error_e=error,
        reserve_a_plus_b=deficit + bridge,
        log_n_minus_last_prime=log_n - last_prime,
        max_exponent=max_exponent,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extend the C=2 Rosser shortfall trend on sampled CA support points."
    )
    parser.add_argument("--max-prime", type=int, default=250_000_000)
    parser.add_argument(
        "--epsilons",
        default=(
            "0.00000002,0.00000001,0.000000005,0.000000002,0.000000001,"
            "0.0000000007,0.0000000005,0.0000000003,0.0000000002"
        ),
    )
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    args = parser.parse_args()

    primes = primes_up_to(args.max_prime)
    prefix_log_beta, prefix_inv_p, prefix_c2, prefix_s_minus1 = build_prefixes(primes)

    points = []
    skipped_epsilons = []
    for epsilon in parse_float_list(args.epsilons):
        point = evaluate_ca_point(primes, prefix_log_beta, epsilon)
        if point is None:
            skipped_epsilons.append(epsilon)
        else:
            points.append(point)
    points = sorted({point.last_prime: point for point in points}.values(), key=lambda item: item.last_prime)

    if len(points) < 2:
        raise ValueError("Need at least two non-truncated CA points.")

    def prime_index(value: int) -> int:
        return bisect.bisect_right(primes, value)

    def a_remainder(value: int) -> float:
        index = prime_index(value)
        return prefix_inv_p[index] - math.log(math.log(value))

    rows: list[C2ShortfallTrendRow] = []
    upper = points[0].beta_error_e
    for previous, current in zip(points, points[1:]):
        y = previous.last_prime
        x = current.last_prime
        y_index = prime_index(y)
        x_index = prime_index(x)
        nu = x_index - y_index
        h_log = math.log(math.log(x) / math.log(y))
        mu = h_log / nu
        exp_mu = math.exp(mu)
        exp_minus_mu = math.exp(-mu)
        s_minus1 = prefix_s_minus1[x_index] - prefix_s_minus1[y_index]
        c2 = prefix_c2[x_index] - prefix_c2[y_index]
        moment_1 = exp_minus_mu * (nu + s_minus1) - nu
        previous_upper = upper
        threshold = current.reserve_a_plus_b - previous_upper
        upper = previous_upper + moment_1
        d_term = nu * (1.0 - (1.0 + mu) * exp_minus_mu)
        required_a_upper = a_remainder(y) - c2 + exp_mu * (threshold + d_term)
        actual_a_x = a_remainder(x)
        scale = math.sqrt(x) * math.log(x)
        log_x = math.log(x)
        required_constant = (required_a_upper - MEISSEL_MERTENS_PRIME_CONSTANT) * scale
        actual_constant = (actual_a_x - MEISSEL_MERTENS_PRIME_CONSTANT) * scale
        finite_c2_constant = 2.0
        c2_shortfall = 2.0 - required_constant
        additive_shortfall = c2_shortfall / scale
        dusart_2018_constant = 0.2 * math.sqrt(x) / (log_x**2)
        axler_2018_constant = None
        axler_2018_margin = None
        if x >= 46_909_074:
            axler_2018_constant = math.sqrt(x) * (
                1.0 / (20.0 * log_x**2) + 3.0 / (16.0 * log_x**3)
            )
            axler_2018_margin = required_constant - axler_2018_constant
        rows.append(
            C2ShortfallTrendRow(
                y_prime=y,
                x_prime=x,
                block_prime_count=nu,
                previous_upper_e_y=previous_upper,
                reserve_x=current.reserve_a_plus_b,
                moment_1=moment_1,
                m1_step_threshold=threshold,
                required_constant=required_constant,
                actual_constant=actual_constant,
                finite_c2_constant=finite_c2_constant,
                c2_shortfall=c2_shortfall,
                additive_shortfall=additive_shortfall,
                dusart_2018_constant=dusart_2018_constant,
                axler_2018_constant=axler_2018_constant,
                axler_2018_margin=axler_2018_margin,
                actual_margin=required_a_upper - actual_a_x,
                cumulative_upper_e_x=upper,
                beta_error_e_x=current.beta_error_e,
                upper_minus_actual_e=upper - current.beta_error_e,
            )
        )

    payload = {
        "notes": [
            "This is a numerical trend audit, not an analytic proof.",
            "finite_c2_constant is the Rosser-Schoenfeld finite-window scale B+2/(sqrt(x)log x), used here as a diagnostic comparator.",
            "c2_shortfall = 2 - C_required. Positive means the C=2 comparator is still too loose by additive_shortfall.",
            "additive_shortfall = (2-C_required)/(sqrt(x)log x), the raw A(x)-scale gap.",
            "dusart_2018_constant converts |A(x)-B| <= 0.2/log^3(x), valid for x >= 2,278,383, into the common C/(sqrt(x)log x) scale.",
            "axler_2018_constant converts Axler's upper bound 1/(20log^3 x)+3/(16log^4 x), valid for x >= 46,909,074, into the common C scale.",
            "Positive axler_2018_margin means Axler's bound is strong enough for that CA endpoint.",
        ],
        "args": {
            "max_prime": args.max_prime,
            "epsilons": parse_float_list(args.epsilons),
            "skipped_epsilons": skipped_epsilons,
        },
        "ca_points": [asdict(point) for point in points],
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    with args.csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))

    print(f"primes up to {args.max_prime}: {len(primes)}")
    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    if skipped_epsilons:
        print(f"skipped truncated epsilons: {skipped_epsilons}")
    print()
    print("C=2 shortfall trend:")
    for row in rows:
        print(
            f"  {row.y_prime:>9}->{row.x_prime:<9} "
            f"C_req={row.required_constant:.6f} "
            f"2-C={row.c2_shortfall:+.6f} "
            f"Axler={row.axler_2018_constant if row.axler_2018_constant is not None else float('nan'):.6f} "
            f"add={row.additive_shortfall:.3e} "
            f"actual_margin={row.actual_margin:+.3e} "
            f"U-E={row.upper_minus_actual_e:+.3e}"
        )


if __name__ == "__main__":
    main()
