from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_THETA_CSV = ROOT / "ca_theta_cancellation_rows.csv"
DEFAULT_JSON = ROOT / "w_mvdc_certified_envelope_report.json"
DEFAULT_CSV = ROOT / "w_mvdc_certified_envelope_rows.csv"


@dataclass
class CertifiedBlockRow:
    y_prime: int
    x_prime: int
    order: int
    endpoint_mode: str
    exact_nu: int
    exact_s1: float
    exact_s2: float
    nu_lower: float
    nu_upper: float
    s1_lower: float
    s1_upper: float
    s2_lower: float
    s2_upper: float
    lambda_upper: float
    p_upper: float
    remainder_upper: float
    delta_rho_lower: float
    w_block_upper: float
    c_step_certified: float | None
    exact_w_block: float
    exact_c_step: float
    certificate_margin: float


def read_theta_rows(path: Path) -> dict[int, dict[str, float]]:
    rows: dict[int, dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            x = int(raw["last_prime"])
            rows[x] = {key: float(value) for key, value in raw.items() if key != "last_prime"}
    return rows


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


def prime_indices(primes: list[int]) -> dict[int, int]:
    return {prime: index for index, prime in enumerate(primes)}


def f_scale(x: int) -> float:
    return 1.0 / (math.sqrt(x) * math.log(x))


def pi_lower(t: float) -> float:
    return t / (math.log(t) - 1.0)


def pi_upper(t: float) -> float:
    return t / (math.log(t) - 1.1)


def integral_t_over_tminus1_power(y: float, x: float, s: int) -> float:
    # Integral of t/(t-1)^(s+1) dt from y to x.
    uy = y - 1.0
    ux = x - 1.0
    if s == 1:
        return (math.log(ux) - 1.0 / ux) - (math.log(uy) - 1.0 / uy)
    primitive_x = (ux ** (1 - s)) / (1 - s) - (ux ** (-s)) / s
    primitive_y = (uy ** (1 - s)) / (1 - s) - (uy ** (-s)) / s
    return primitive_x - primitive_y


def s_bounds(y: int, x: int, s: int) -> tuple[float, float]:
    if s == 0:
        lower = pi_lower(x) - pi_upper(y)
        upper = pi_upper(x) - pi_lower(y)
        return max(0.0, lower), max(0.0, upper)

    f_x = 1.0 / ((x - 1.0) ** s)
    f_y = 1.0 / ((y - 1.0) ** s)
    integral = integral_t_over_tminus1_power(float(y), float(x), s)

    lower = f_x * pi_lower(x) - f_y * pi_upper(y) + (s / (math.log(x) - 1.0)) * integral
    upper = f_x * pi_upper(x) - f_y * pi_lower(y) + (s / (math.log(y) - 1.1)) * integral
    return max(0.0, lower), max(0.0, upper)


def exact_block_values(primes: list[int], indices: dict[int, int], y: int, x: int) -> tuple[int, float, float, float, float]:
    y_index = indices[y]
    x_index = indices[x]
    block = primes[y_index + 1 : x_index + 1]
    main = math.log(math.log(x) / math.log(y))
    q_value = sum(-math.log1p(-1.0 / prime) for prime in block) - main
    s1 = sum(1.0 / (prime - 1.0) for prime in block)
    s2 = sum(1.0 / ((prime - 1.0) ** 2) for prime in block)
    return len(block), main, q_value, s1, s2


def reduced_moment_interval(
    r: int,
    mu: float,
    lower_s: dict[int, float],
    upper_s: dict[int, float],
) -> tuple[float, float]:
    base = math.exp(-mu) - 1.0
    coeff = math.exp(-mu)
    low = 0.0
    high = 0.0
    for s in range(0, r + 1):
        a = math.comb(r, s) * (base ** (r - s)) * (coeff**s)
        if a >= 0.0:
            low += a * lower_s[s]
            high += a * upper_s[s]
        else:
            low += a * upper_s[s]
            high += a * lower_s[s]
    return low, high


def p_upper_from_moment_intervals(moment_intervals: dict[int, tuple[float, float]], order: int) -> float:
    total = 0.0
    for r in range(1, order + 1):
        coeff = ((-1.0) ** (r - 1)) / r
        low, high = moment_intervals[r]
        total += coeff * (high if coeff >= 0.0 else low)
    return total


def lambda_upper(y: int, x: int, mu: float) -> float:
    base = math.exp(-mu) - 1.0
    coeff = math.exp(-mu)
    # For p>Y, 1/(p-1) < 1/(Y-1), and for p<=x it is >= 1/(x-1).
    upper_endpoint = base + coeff / (y - 1.0)
    lower_endpoint = base + coeff / (x - 1.0)
    return max(abs(upper_endpoint), abs(lower_endpoint))


def delta_rho_lower_bound(
    y: int,
    x: int,
    theta_rows: dict[int, dict[str, float]],
    endpoint_mode: str,
    theta_c: float,
    theta_power: int,
) -> float:
    if endpoint_mode == "exact":
        return theta_rows[x]["endpoint_term"] - theta_rows[y]["endpoint_term"]
    if endpoint_mode != "theta-bound":
        raise ValueError(f"Unknown endpoint mode: {endpoint_mode}")
    exponent = theta_power + 1
    rho_x_lower = -theta_c / (math.log(x) ** exponent)
    rho_y_upper = theta_c / (math.log(y) ** exponent)
    return rho_x_lower - rho_y_upper


def evaluate_pair(
    primes: list[int],
    indices: dict[int, int],
    theta_rows: dict[int, dict[str, float]],
    y: int,
    x: int,
    order: int,
    endpoint_mode: str,
    theta_c: float,
    theta_power: int,
) -> CertifiedBlockRow:
    exact_nu, main, q_value, exact_s1, exact_s2 = exact_block_values(primes, indices, y, x)
    mu = main / exact_nu

    lower_s: dict[int, float] = {}
    upper_s: dict[int, float] = {}
    for s in range(0, order + 1):
        lower_s[s], upper_s[s] = s_bounds(y, x, s)

    moment_intervals = {
        r: reduced_moment_interval(r, mu, lower_s, upper_s)
        for r in range(1, order + 1)
    }
    lam = lambda_upper(y, x, mu)
    remainder = upper_s[0] * (lam ** (order + 1)) / ((order + 1) * (1.0 - lam))
    p_upper = p_upper_from_moment_intervals(moment_intervals, order)
    delta_lower = delta_rho_lower_bound(y, x, theta_rows, endpoint_mode, theta_c, theta_power)
    w_upper = p_upper + remainder - delta_lower

    delta_exact = theta_rows[x]["endpoint_term"] - theta_rows[y]["endpoint_term"]
    w_exact = q_value - delta_exact
    f_gap = f_scale(y) - f_scale(x)
    c_certified = (-w_upper / f_gap) if w_upper < 0.0 else None
    c_exact = -w_exact / f_gap

    return CertifiedBlockRow(
        y_prime=y,
        x_prime=x,
        order=order,
        endpoint_mode=endpoint_mode,
        exact_nu=exact_nu,
        exact_s1=exact_s1,
        exact_s2=exact_s2,
        nu_lower=lower_s[0],
        nu_upper=upper_s[0],
        s1_lower=lower_s.get(1, 0.0),
        s1_upper=upper_s.get(1, 0.0),
        s2_lower=lower_s.get(2, 0.0),
        s2_upper=upper_s.get(2, 0.0),
        lambda_upper=lam,
        p_upper=p_upper,
        remainder_upper=remainder,
        delta_rho_lower=delta_lower,
        w_block_upper=w_upper,
        c_step_certified=c_certified,
        exact_w_block=w_exact,
        exact_c_step=c_exact,
        certificate_margin=w_upper - w_exact,
    )


def write_csv(path: Path, rows: list[CertifiedBlockRow]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="First certified-envelope audit for corrected MVDC W-blocks.")
    parser.add_argument("--theta-csv", type=Path, default=DEFAULT_THETA_CSV)
    parser.add_argument("--max-prime", type=int, default=120_000_000)
    parser.add_argument("--min-prime", type=int, default=3_329_267)
    parser.add_argument("--order", type=int, default=6)
    parser.add_argument("--endpoint-mode", choices=["exact", "theta-bound"], default="exact")
    parser.add_argument("--theta-c", type=float, default=0.006788)
    parser.add_argument("--theta-power", type=int, default=1)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    theta_rows = read_theta_rows(args.theta_csv)
    targets = [prime for prime in sorted(theta_rows) if args.min_prime <= prime <= args.max_prime]
    if len(targets) < 2:
        raise ValueError("Need at least two theta rows in range")

    primes = primes_up_to(max(args.max_prime, targets[-1]))
    indices = prime_indices(primes)
    rows = [
        evaluate_pair(
            primes=primes,
            indices=indices,
            theta_rows=theta_rows,
            y=y,
            x=x,
            order=args.order,
            endpoint_mode=args.endpoint_mode,
            theta_c=args.theta_c,
            theta_power=args.theta_power,
        )
        for y, x in zip(targets, targets[1:])
    ]

    payload = {
        "notes": [
            "This is a first certified-envelope audit for corrected MVDC W-blocks.",
            "Pi bounds use Dusart-style inputs pi(t)>t/(log t-1), pi(t)<t/(log t-1.1) in the sampled range.",
            "endpoint-mode=exact isolates the S_s certification loss; endpoint-mode=theta-bound also applies a symmetric theta envelope.",
            "The theta-bound mode is expected to be very coarse and is a diagnostic of endpoint-control difficulty.",
            "The present audit fixes mu from the exact finite block count; an infinite-range proof would also need either a certified nu mesh or a mu-interval version.",
            "The current loss is dominated by independent S0/S1 interval bounds, which destroy the first-moment cancellation.",
        ],
        "args": {
            "theta_csv": str(args.theta_csv),
            "max_prime": args.max_prime,
            "min_prime": args.min_prime,
            "order": args.order,
            "endpoint_mode": args.endpoint_mode,
            "theta_c": args.theta_c,
            "theta_power": args.theta_power,
        },
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    write_csv(args.csv, rows)

    print(f"primes up to {max(args.max_prime, targets[-1])}: {len(primes)}")
    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print()
    print(f"Certified corrected W MVDC envelope ({args.endpoint_mode} endpoint):")
    for row in rows:
        cert = row.c_step_certified if row.c_step_certified is not None else float("nan")
        print(
            f"{row.y_prime:>9}->{row.x_prime:<9} "
            f"C_exact={row.exact_c_step:.6f} C_cert={cert:.6f} "
            f"margin={row.certificate_margin:+.3e} "
            f"nu=[{row.nu_lower:.1f},{row.nu_upper:.1f}] "
            f"S1=[{row.s1_lower:.6e},{row.s1_upper:.6e}]"
        )


if __name__ == "__main__":
    main()
