from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable

from test_w_mvdc_certified_envelope import DEFAULT_THETA_CSV, f_scale, primes_up_to


ROOT = Path(__file__).resolve().parent
DEFAULT_JSON = ROOT / "mvdc_signed_chebyshev_center_report.json"
DEFAULT_CSV = ROOT / "mvdc_signed_chebyshev_center_rows.csv"


@dataclass
class SignedChebyshevCenterRow:
    y_prime: int
    x_prime: int
    prime_count: int
    chebyshev_integral: float
    sqrt_weight_integral: float
    mvdc_center_k: float
    required_k_for_c_190: float
    required_k_for_c_200: float
    required_k_for_c_205: float
    center_margin_c_200: float
    weighted_variance: float
    weighted_std: float
    weighted_skewness: float
    weighted_kurtosis: float
    min_sample_k: float
    max_sample_k: float
    mean_sample_k: float
    s2_block: float
    f_gap: float


def read_targets(path: Path, min_prime: int, max_prime: int) -> list[int]:
    targets: list[int] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            prime = int(raw["last_prime"])
            if min_prime <= prime <= max_prime:
                targets.append(prime)
    return sorted(targets)


def adaptive_simpson(
    fn: Callable[[float], float],
    left: float,
    right: float,
    eps: float = 1e-13,
    max_depth: int = 30,
) -> float:
    def simpson(a: float, b: float) -> float:
        c = 0.5 * (a + b)
        return (b - a) * (fn(a) + 4.0 * fn(c) + fn(b)) / 6.0

    def recurse(a: float, b: float, whole: float, depth: int) -> float:
        c = 0.5 * (a + b)
        left_part = simpson(a, c)
        right_part = simpson(c, b)
        total = left_part + right_part
        if depth <= 0 or abs(total - whole) <= 15.0 * eps:
            return total + (total - whole) / 15.0
        return recurse(a, c, left_part, depth - 1) + recurse(c, b, right_part, depth - 1)

    whole_interval = simpson(left, right)
    return recurse(left, right, whole_interval, max_depth)


def sqrt_weight_integral(y: int, x: int) -> float:
    left = math.log(y)
    right = math.log(x)

    def integrand(u: float) -> float:
        return (u + 1.0) * math.exp(-0.5 * u) / (u * u)

    return adaptive_simpson(integrand, left, right)


def g_weight(t: float) -> float:
    return 1.0 / (t * math.log(t))


def t_omega_antiderivative(t: float) -> float:
    return math.log(math.log(t)) - 1.0 / math.log(t)


def chebyshev_interval_integral(theta_value: float, left: float, right: float) -> float:
    return theta_value * (g_weight(left) - g_weight(right)) - (
        t_omega_antiderivative(right) - t_omega_antiderivative(left)
    )


def sqrt_weight_density(t: float) -> float:
    log_t = math.log(t)
    return (log_t + 1.0) / (t ** 1.5 * log_t * log_t)


def s2_block(primes: list[int], left_index: int, right_index: int) -> float:
    return sum(-math.log1p(-1.0 / prime) - 1.0 / prime for prime in primes[left_index + 1 : right_index + 1])


def required_k(c_value: float, f_gap: float, s2_value: float, sqrt_weight: float) -> float:
    return (c_value * f_gap + s2_value) / sqrt_weight


def evaluate_block(
    primes: list[int],
    theta_prefix: list[float],
    y_index: int,
    x_index: int,
    y: int,
    x: int,
) -> SignedChebyshevCenterRow:
    cheb_integral = 0.0
    samples: list[tuple[float, float]] = []
    current_left = float(y)
    theta_value = theta_prefix[y_index]

    for index in range(y_index + 1, x_index + 1):
        right = float(primes[index])
        if right > current_left:
            cheb_integral += chebyshev_interval_integral(theta_value, current_left, right)
            mid = 0.5 * (current_left + right)
            width = right - current_left
            weight = sqrt_weight_density(mid) * width
            k_value = -(theta_value - mid) / math.sqrt(mid)
            samples.append((weight, k_value))
        theta_value = theta_prefix[index]
        current_left = right

    sqrt_weight = sqrt_weight_integral(y, x)
    center_k = -cheb_integral / sqrt_weight
    total_weight_sample = sum(weight for weight, _value in samples)
    mean_sample = sum(weight * value for weight, value in samples) / total_weight_sample
    centered = [(weight, value - center_k) for weight, value in samples]
    variance = sum(weight * delta * delta for weight, delta in centered) / total_weight_sample
    std = math.sqrt(max(0.0, variance))
    if std > 0.0:
        skewness = sum(weight * (delta / std) ** 3 for weight, delta in centered) / total_weight_sample
        kurtosis = sum(weight * (delta / std) ** 4 for weight, delta in centered) / total_weight_sample
    else:
        skewness = 0.0
        kurtosis = 0.0
    min_sample = min(value for _weight, value in samples)
    max_sample = max(value for _weight, value in samples)
    s2_value = s2_block(primes, y_index, x_index)
    f_gap = f_scale(y) - f_scale(x)
    k190 = required_k(1.90, f_gap, s2_value, sqrt_weight)
    k200 = required_k(2.00, f_gap, s2_value, sqrt_weight)
    k205 = required_k(2.05, f_gap, s2_value, sqrt_weight)
    return SignedChebyshevCenterRow(
        y_prime=y,
        x_prime=x,
        prime_count=x_index - y_index,
        chebyshev_integral=cheb_integral,
        sqrt_weight_integral=sqrt_weight,
        mvdc_center_k=center_k,
        required_k_for_c_190=k190,
        required_k_for_c_200=k200,
        required_k_for_c_205=k205,
        center_margin_c_200=center_k - k200,
        weighted_variance=variance,
        weighted_std=std,
        weighted_skewness=skewness,
        weighted_kurtosis=kurtosis,
        min_sample_k=min_sample,
        max_sample_k=max_sample,
        mean_sample_k=mean_sample,
        s2_block=s2_value,
        f_gap=f_gap,
    )


def prime_indices(primes: list[int]) -> dict[int, int]:
    return {prime: index for index, prime in enumerate(primes)}


def theta_prefixes(primes: list[int]) -> list[float]:
    theta = 0.0
    values: list[float] = []
    for prime in primes:
        theta += math.log(prime)
        values.append(theta)
    return values


def write_csv(path: Path, rows: list[SignedChebyshevCenterRow]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="MVDC-centred profile of the signed Chebyshev block needed for W.")
    parser.add_argument("--theta-csv", type=Path, default=DEFAULT_THETA_CSV)
    parser.add_argument("--min-prime", type=int, default=3_329_267)
    parser.add_argument("--max-prime", type=int, default=56_048_351)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    targets = read_targets(args.theta_csv, args.min_prime, args.max_prime)
    if len(targets) < 2:
        raise ValueError("Need at least two target primes")
    primes = primes_up_to(max(args.max_prime, targets[-1]))
    indices = prime_indices(primes)
    theta_values = theta_prefixes(primes)
    rows = [
        evaluate_block(primes, theta_values, indices[y], indices[x], y, x)
        for y, x in zip(targets, targets[1:])
    ]
    payload = {
        "notes": [
            "This is the MVDC centre profile for the signed Chebyshev block.",
            "The centre K is defined by integral(theta(t)-t) omega(t) dt = -K integral sqrt(t) omega(t) dt.",
            "This does not use Dusart/Axler pi envelopes; it studies the actual centred object that the W-ledger needs.",
            "Moment columns are diagnostic midpoint quadrature for the centred K-profile; the first moment is zero by construction of K.",
        ],
        "args": {
            "theta_csv": str(args.theta_csv),
            "min_prime": args.min_prime,
            "max_prime": args.max_prime,
        },
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    write_csv(args.csv, rows)

    print(f"primes up to {max(args.max_prime, targets[-1])}: {len(primes)}")
    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print()
    print("MVDC signed Chebyshev centres:")
    for row in rows:
        print(
            f"{row.y_prime:>9}->{row.x_prime:<9} "
            f"K={row.mvdc_center_k:.6f} Kreq2={row.required_k_for_c_200:.6f} "
            f"margin={row.center_margin_c_200:+.6f} "
            f"std={row.weighted_std:.6f} skew={row.weighted_skewness:+.3f} "
            f"range=[{row.min_sample_k:.3f},{row.max_sample_k:.3f}]"
        )


if __name__ == "__main__":
    main()
