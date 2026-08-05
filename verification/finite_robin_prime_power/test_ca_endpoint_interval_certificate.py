from __future__ import annotations

import argparse
import decimal
import hashlib
import json
import math
import sys
import time
from dataclasses import asdict, dataclass
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN
from pathlib import Path
from typing import Iterator


ROOT = Path(__file__).resolve().parent
DEFAULT_REPORT = ROOT / "ca_endpoint_interval_certificate_report.json"


@dataclass(frozen=True)
class EndpointSpec:
    epsilon_text: str
    support_prime: int


@dataclass(frozen=True)
class BlockSpec:
    y_prime: int
    x_prime: int
    k_upper_text: str


ENDPOINTS = (
    EndpointSpec("0.00000002", 3_329_267),
    EndpointSpec("0.00000001", 6_382_007),
    EndpointSpec("0.000000005", 12_253_883),
    EndpointSpec("0.000000002", 29_093_377),
    EndpointSpec("0.000000001", 56_048_351),
)

BLOCKS = (
    BlockSpec(3_329_267, 6_382_007, "0.15"),
    BlockSpec(6_382_007, 12_253_883, "0.14"),
    BlockSpec(12_253_883, 29_093_377, "0.40"),
    BlockSpec(29_093_377, 56_048_351, "0.11"),
)


Interval = tuple[Decimal, Decimal]


class DecimalIntervals:
    """Small outward-rounded interval layer over CPython Decimal."""

    def __init__(self, precision: int) -> None:
        if precision < 30:
            raise ValueError("precision must be at least 30 decimal digits")
        self.precision = precision
        self.nearest = Context(prec=precision, rounding=ROUND_HALF_EVEN)
        self.floor = Context(prec=precision, rounding=ROUND_FLOOR)
        self.ceiling = Context(prec=precision, rounding=ROUND_CEILING)
        self.padding = Context(prec=precision + 8, rounding=ROUND_HALF_EVEN)
        self._ulp_cache: dict[int, Decimal] = {}

    def point(self, value: int | str | Decimal) -> Interval:
        number = value if isinstance(value, Decimal) else Decimal(value)
        return number, number

    def _ulp(self, value: Decimal) -> Decimal:
        exponent = value.adjusted() - self.precision + 1
        if exponent not in self._ulp_cache:
            self._ulp_cache[exponent] = Decimal(f"1e{exponent}")
        return self._ulp_cache[exponent]

    def _padded_unary(self, operation: str, value: Decimal) -> Interval:
        if operation == "ln":
            midpoint = self.nearest.ln(value)
        elif operation == "exp":
            midpoint = self.nearest.exp(value)
        elif operation == "sqrt":
            midpoint = self.nearest.sqrt(value)
        else:
            raise ValueError(f"unsupported unary operation {operation!r}")
        ulp = self._ulp(midpoint)
        return (
            self.padding.subtract(midpoint, ulp),
            self.padding.add(midpoint, ulp),
        )

    def ln_integer(self, value: int) -> Interval:
        if value <= 0:
            raise ValueError("logarithm input must be positive")
        return self._padded_unary("ln", Decimal(value))

    def ln(self, interval: Interval) -> Interval:
        lower, upper = interval
        if lower <= 0:
            raise ValueError("logarithm interval must be positive")
        lower_result = self._padded_unary("ln", lower)[0]
        upper_result = self._padded_unary("ln", upper)[1]
        return lower_result, upper_result

    def exp(self, interval: Interval) -> Interval:
        lower_result = self._padded_unary("exp", interval[0])[0]
        upper_result = self._padded_unary("exp", interval[1])[1]
        return lower_result, upper_result

    def sqrt_integer(self, value: int) -> Interval:
        if value < 0:
            raise ValueError("square-root input must be nonnegative")
        return self._padded_unary("sqrt", Decimal(value))

    def add(self, left: Interval, right: Interval) -> Interval:
        return (
            self.floor.add(left[0], right[0]),
            self.ceiling.add(left[1], right[1]),
        )

    def sub(self, left: Interval, right: Interval) -> Interval:
        return (
            self.floor.subtract(left[0], right[1]),
            self.ceiling.subtract(left[1], right[0]),
        )

    def mul(self, left: Interval, right: Interval) -> Interval:
        lower_values = [
            self.floor.multiply(a, b) for a in left for b in right
        ]
        upper_values = [
            self.ceiling.multiply(a, b) for a in left for b in right
        ]
        return min(lower_values), max(upper_values)

    def div(self, numerator: Interval, denominator: Interval) -> Interval:
        if denominator[0] <= 0 <= denominator[1]:
            raise ZeroDivisionError("division interval contains zero")
        lower_values = [
            self.floor.divide(a, b) for a in numerator for b in denominator
        ]
        upper_values = [
            self.ceiling.divide(a, b) for a in numerator for b in denominator
        ]
        return min(lower_values), max(upper_values)

    def rational(self, numerator: int, denominator: int) -> Interval:
        if denominator <= 0:
            raise ValueError("rational denominator must be positive")
        n_value = Decimal(numerator)
        d_value = Decimal(denominator)
        return (
            self.floor.divide(n_value, d_value),
            self.ceiling.divide(n_value, d_value),
        )

    def scale_integer(self, interval: Interval, multiplier: int) -> Interval:
        return self.mul(interval, self.point(multiplier))


@dataclass
class PrefixValues:
    theta: Interval
    squarefree_log_sigma: Interval
    reciprocal_prime_sum: Interval
    log_prime: Interval


@dataclass
class ProfileAudit:
    epsilon_text: str
    support_prime: int
    next_prime: int
    support_lower: Interval
    support_upper: Interval
    representative_inside_support: bool
    layer_endpoints: list[int]
    layer_checks: list[dict[str, object]]
    no_further_layer_check: dict[str, object]
    unique_on_open_support_interval: bool


def interval_row(interval: Interval) -> dict[str, str]:
    return {"lower": str(interval[0]), "upper": str(interval[1])}


def odd_sieve(limit: int) -> bytearray:
    if limit < 2:
        raise ValueError("sieve limit must be at least 2")
    flags = bytearray(b"\x01") * (limit // 2 + 1)
    flags[0] = 0
    for prime in range(3, math.isqrt(limit) + 1, 2):
        if flags[prime // 2]:
            start = prime * prime // 2
            count = (len(flags) - start - 1) // prime + 1
            flags[start::prime] = b"\x00" * count
    return flags


def is_prime(flags: bytearray, value: int) -> bool:
    if value == 2:
        return True
    return value >= 3 and value % 2 == 1 and bool(flags[value // 2])


def primes_from_sieve(flags: bytearray, limit: int) -> Iterator[int]:
    yield 2
    for value in range(3, limit + 1, 2):
        if flags[value // 2]:
            yield value


def next_prime(flags: bytearray, value: int, limit: int) -> int:
    candidate = value + 1 if value == 2 else value + 2
    if candidate % 2 == 0:
        candidate += 1
    while candidate <= limit:
        if is_prime(flags, candidate):
            return candidate
        candidate += 2
    raise ValueError(f"no next prime after {value} within sieve limit {limit}")


def transition_tau(
    arithmetic: DecimalIntervals, prime: int, layer: int
) -> Interval:
    """Enclose the CA gain threshold for adding the `layer`-th prime factor."""

    # The local gain is positive precisely when epsilon < tau_layer(p), where
    #
    #   tau_layer(p) = log((1-p^(-layer-1))/(1-p^(-layer))) / log(p).
    #
    # The ratio below is formed as an exact rational before taking its log.
    numerator = prime ** (layer + 1) - 1
    denominator = prime * (prime**layer - 1)
    ratio = arithmetic.rational(numerator, denominator)
    return arithmetic.div(arithmetic.ln(ratio), arithmetic.ln_integer(prime))


def transition_tau_float(prime: int, layer: int) -> float:
    inverse_power = prime ** (-layer)
    numerator = math.log1p(-inverse_power / prime) - math.log1p(-inverse_power)
    return numerator / math.log(prime)


def last_transition_above(
    primes: list[int], layer: int, epsilon: float
) -> int:
    low = 0
    high = len(primes)
    while low < high:
        middle = (low + high) // 2
        if transition_tau_float(primes[middle], layer) > epsilon:
            low = middle + 1
        else:
            high = middle
    return low - 1


def audit_profile(
    arithmetic: DecimalIntervals,
    spec: EndpointSpec,
    next_support_prime: int,
    search_primes: list[int],
) -> ProfileAudit:
    epsilon = Decimal(spec.epsilon_text)
    support_lower = transition_tau(arithmetic, next_support_prime, 1)
    support_upper = transition_tau(arithmetic, spec.support_prime, 1)
    representative_inside = support_lower[1] < epsilon < support_upper[0]

    layer_endpoints = [spec.support_prime]
    layer_checks: list[dict[str, object]] = [
        {
            "layer": 1,
            "last_prime": spec.support_prime,
            "next_prime": next_support_prime,
            "tau_at_last": interval_row(support_upper),
            "tau_at_next": interval_row(support_lower),
            "representative_check": representative_inside,
            "constant_across_open_support_interval": True,
        }
    ]

    layer = 2
    all_constant = representative_inside
    no_further: dict[str, object] | None = None
    while True:
        tau_at_two = transition_tau(arithmetic, 2, layer)
        if tau_at_two[1] < support_lower[0]:
            no_further = {
                "first_absent_layer": layer,
                "tau_at_two": interval_row(tau_at_two),
                "below_support_lower": True,
            }
            break

        index = last_transition_above(search_primes, layer, float(epsilon))
        if index < 0 or index + 1 >= len(search_primes):
            raise ValueError(
                f"profile search range is insufficient at layer {layer} for "
                f"support {spec.support_prime}"
            )
        last_prime = search_primes[index]
        following_prime = search_primes[index + 1]
        tau_at_last = transition_tau(arithmetic, last_prime, layer)
        tau_at_next = transition_tau(arithmetic, following_prime, layer)
        representative_check = tau_at_last[0] > epsilon > tau_at_next[1]
        constant_check = (
            tau_at_last[0] > support_upper[1]
            and tau_at_next[1] < support_lower[0]
        )
        all_constant = all_constant and representative_check and constant_check
        layer_endpoints.append(last_prime)
        layer_checks.append(
            {
                "layer": layer,
                "last_prime": last_prime,
                "next_prime": following_prime,
                "tau_at_last": interval_row(tau_at_last),
                "tau_at_next": interval_row(tau_at_next),
                "representative_check": representative_check,
                "constant_across_open_support_interval": constant_check,
            }
        )
        layer += 1

    assert no_further is not None
    return ProfileAudit(
        epsilon_text=spec.epsilon_text,
        support_prime=spec.support_prime,
        next_prime=next_support_prime,
        support_lower=support_lower,
        support_upper=support_upper,
        representative_inside_support=representative_inside,
        layer_endpoints=layer_endpoints,
        layer_checks=layer_checks,
        no_further_layer_check=no_further,
        unique_on_open_support_interval=all_constant,
    )


def harmonic_gamma_interval(
    arithmetic: DecimalIntervals, term_count: int
) -> Interval:
    """Self-certify gamma with elementary harmonic-number inequalities."""

    harmonic: Interval = arithmetic.point(0)
    for denominator in range(1, term_count + 1):
        harmonic = arithmetic.add(
            harmonic, arithmetic.rational(1, denominator)
        )
    log_n = arithmetic.ln_integer(term_count)
    # 1/(2n+1) < H_n-log(n)-gamma < 1/(2n).
    lower = arithmetic.sub(
        arithmetic.sub(harmonic, log_n),
        arithmetic.rational(1, 2 * term_count),
    )[0]
    upper = arithmetic.sub(
        arithmetic.sub(harmonic, log_n),
        arithmetic.rational(1, 2 * term_count + 1),
    )[1]
    return lower, upper


def squarefree_sigma_log_factor(
    arithmetic: DecimalIntervals, prime: int, series_cutoff: int
) -> Interval:
    """Enclose log(1+1/p); use an alternating rational tail for large p."""

    if prime <= series_cutoff:
        return arithmetic.ln(arithmetic.rational(prime + 1, prime))

    first = arithmetic.rational(1, prime)
    second = arithmetic.rational(1, 2 * prime * prime)
    third = arithmetic.rational(1, 3 * prime * prime * prime)
    # S_2 <= log(1+1/p) <= S_3 for the alternating logarithm series.
    lower = arithmetic.sub(first, second)[0]
    upper = arithmetic.add(
        arithmetic.sub(first, second), third
    )[1]
    return lower, upper


def sigma_profile_correction(
    arithmetic: DecimalIntervals, prime: int, exponent: int
) -> Interval:
    """Enclose log((sigma(p^a)/p^a)/(1+1/p)) for a >= 2."""

    if exponent < 2:
        return arithmetic.point(0)
    numerator = prime ** (exponent + 1) - 1
    denominator = prime ** (exponent - 1) * (prime * prime - 1)
    return arithmetic.ln(arithmetic.rational(numerator, denominator))


def scan_prime_prefixes(
    arithmetic: DecimalIntervals,
    sieve: bytearray,
    max_prime: int,
    endpoint_set: set[int],
    theta_checkpoint_set: set[int],
    series_cutoff: int,
    progress_every: int,
) -> tuple[
    dict[int, PrefixValues],
    dict[int, Interval],
    int,
    str,
]:
    theta: Interval = arithmetic.point(0)
    squarefree_log_sigma: Interval = arithmetic.point(0)
    reciprocal_sum: Interval = arithmetic.point(0)
    endpoint_values: dict[int, PrefixValues] = {}
    theta_checkpoints: dict[int, Interval] = {}
    prime_count = 0
    digest = hashlib.sha256()
    started = time.perf_counter()

    for prime in primes_from_sieve(sieve, max_prime):
        log_prime = arithmetic.ln_integer(prime)
        theta = arithmetic.add(theta, log_prime)
        squarefree_log_sigma = arithmetic.add(
            squarefree_log_sigma,
            squarefree_sigma_log_factor(arithmetic, prime, series_cutoff),
        )
        reciprocal_sum = arithmetic.add(
            reciprocal_sum, arithmetic.rational(1, prime)
        )
        prime_count += 1
        digest.update(prime.to_bytes(8, byteorder="little", signed=False))

        if prime in theta_checkpoint_set:
            theta_checkpoints[prime] = theta
        if prime in endpoint_set:
            endpoint_values[prime] = PrefixValues(
                theta=theta,
                squarefree_log_sigma=squarefree_log_sigma,
                reciprocal_prime_sum=reciprocal_sum,
                log_prime=log_prime,
            )
        if progress_every and prime_count % progress_every == 0:
            elapsed = time.perf_counter() - started
            print(
                f"processed {prime_count:,} primes through {prime:,} "
                f"in {elapsed:.1f}s",
                flush=True,
            )

    return endpoint_values, theta_checkpoints, prime_count, digest.hexdigest()


def exponent_from_layers(prime: int, layer_endpoints: list[int]) -> int:
    return sum(1 for endpoint in layer_endpoints if prime <= endpoint)


def endpoint_certificate(
    arithmetic: DecimalIntervals,
    audit: ProfileAudit,
    endpoint_prefix: PrefixValues,
    theta_checkpoints: dict[int, Interval],
    small_primes: list[int],
    gamma: Interval,
) -> tuple[dict[str, object], Interval, Interval]:
    log_n = endpoint_prefix.theta
    for layer_endpoint in audit.layer_endpoints[1:]:
        log_n = arithmetic.add(log_n, theta_checkpoints[layer_endpoint])

    log_sigma = endpoint_prefix.squarefree_log_sigma
    second_layer_endpoint = (
        audit.layer_endpoints[1] if len(audit.layer_endpoints) > 1 else 1
    )
    corrected_prime_count = 0
    max_exponent = 1
    for prime in small_primes:
        if prime > second_layer_endpoint:
            break
        exponent = exponent_from_layers(prime, audit.layer_endpoints)
        if exponent >= 2:
            correction = sigma_profile_correction(arithmetic, prime, exponent)
            log_sigma = arithmetic.add(log_sigma, correction)
            corrected_prime_count += 1
            max_exponent = max(max_exponent, exponent)

    log_log_n = arithmetic.ln(log_n)
    log_log_log_n = arithmetic.ln(log_log_n)
    gap = arithmetic.sub(arithmetic.add(gamma, log_log_log_n), log_sigma)
    scale = arithmetic.mul(
        arithmetic.sqrt_integer(audit.support_prime), endpoint_prefix.log_prime
    )
    scaled_gap = arithmetic.mul(gap, scale)

    exponent_runs = [
        {
            "lower_prime_exclusive": (
                1
                if layer + 1 == len(audit.layer_endpoints)
                else audit.layer_endpoints[layer + 1]
            ),
            "upper_prime_inclusive": endpoint,
            "exponent_exact": layer + 1,
        }
        for layer, endpoint in enumerate(audit.layer_endpoints)
    ]
    row: dict[str, object] = {
        "support_prime": audit.support_prime,
        "chosen_epsilon_exact_decimal": audit.epsilon_text,
        "next_prime": audit.next_prime,
        "support_parameter_interval": {
            "open_lower_tau1_next": interval_row(audit.support_lower),
            "closed_upper_tau1_support": interval_row(audit.support_upper),
        },
        "chosen_epsilon_strictly_inside": audit.representative_inside_support,
        "unique_profile_on_open_support_interval": audit.unique_on_open_support_interval,
        "unique_profile_for_exact_support_on_closed_interval": (
            audit.unique_on_open_support_interval
        ),
        "layer_endpoints": audit.layer_endpoints,
        "max_exponent": max_exponent,
        "primes_with_exponent_at_least_two": corrected_prime_count,
        "exponent_layer_description": exponent_runs,
        "profile_transition_checks": audit.layer_checks,
        "no_further_layer_check": audit.no_further_layer_check,
        "log_n_interval": interval_row(log_n),
        "log_sigma_over_n_interval": interval_row(log_sigma),
        "exact_log_robin_gap_interval": interval_row(gap),
        "scaled_gap_interval": interval_row(scaled_gap),
        "status": "PASS" if gap[0] > 0 else "FAIL",
    }
    return row, gap, endpoint_prefix.theta


def rho_interval(
    arithmetic: DecimalIntervals,
    endpoint: int,
    prefix: PrefixValues,
) -> Interval:
    numerator = arithmetic.sub(prefix.theta, arithmetic.point(endpoint))
    denominator = arithmetic.mul(
        arithmetic.point(endpoint), prefix.log_prime
    )
    return arithmetic.div(numerator, denominator)


def w0_lower_bound(
    arithmetic: DecimalIntervals,
    y_prime: int,
    x_prime: int,
    steps: int,
) -> Decimal:
    """Lower-bound W0 by right rectangles after the substitution z=log(t)."""

    if steps <= 0:
        raise ValueError("W0 rectangle count must be positive")
    log_y = arithmetic.ln_integer(y_prime)
    log_x = arithmetic.ln_integer(x_prime)
    # [log_y.upper, log_x.lower] is a certified inner integration interval.
    inner_width = arithmetic.floor.subtract(log_x[0], log_y[1])
    h = arithmetic.floor.divide(inner_width, Decimal(steps))
    total = Decimal(0)
    for index in range(1, steps + 1):
        displacement = arithmetic.ceiling.multiply(h, Decimal(index))
        z = arithmetic.ceiling.add(log_y[1], displacement)
        minus_half_z = arithmetic.point(
            arithmetic.floor.divide(arithmetic.floor.minus(z), Decimal(2))
        )
        exponential_lower = arithmetic.exp(minus_half_z)[0]
        z_interval = arithmetic.point(z)
        rational_factor = arithmetic.div(
            arithmetic.add(z_interval, arithmetic.point(1)),
            arithmetic.mul(z_interval, z_interval),
        )[0]
        function_lower = arithmetic.floor.multiply(
            exponential_lower, rational_factor
        )
        total = arithmetic.floor.add(total, function_lower)
    return arithmetic.floor.multiply(h, total)


def block_certificate(
    arithmetic: DecimalIntervals,
    block: BlockSpec,
    prefix_values: dict[int, PrefixValues],
    gap_by_endpoint: dict[int, Interval],
    w0_steps: int,
) -> dict[str, object]:
    y_prefix = prefix_values[block.y_prime]
    x_prefix = prefix_values[block.x_prime]
    reciprocal_block = arithmetic.sub(
        x_prefix.reciprocal_prime_sum, y_prefix.reciprocal_prime_sum
    )
    log_log_ratio = arithmetic.sub(
        arithmetic.ln(x_prefix.log_prime), arithmetic.ln(y_prefix.log_prime)
    )
    rho_y = rho_interval(arithmetic, block.y_prime, y_prefix)
    rho_x = rho_interval(arithmetic, block.x_prime, x_prefix)

    # With U=W(Y) and R=W(x)+G(x), the exact K_req numerator is
    #
    #   log(log x/log Y) - sum_{Y<p<=x} 1/p
    #   + rho(x)-rho(Y)-G(x).
    #
    # The second-order product tail cancels algebraically, so no numerical
    # tail approximation is present in this certificate.
    numerator = arithmetic.sub(log_log_ratio, reciprocal_block)
    numerator = arithmetic.add(numerator, arithmetic.sub(rho_x, rho_y))
    numerator = arithmetic.sub(numerator, gap_by_endpoint[block.x_prime])
    w0_lower = w0_lower_bound(
        arithmetic, block.y_prime, block.x_prime, w0_steps
    )
    if w0_lower <= 0:
        raise ValueError("computed W0 lower bound is not positive")
    if numerator[0] <= 0:
        raise ValueError(
            "the one-sided K_req quotient requires a positive numerator interval"
        )
    k_req_upper = arithmetic.ceiling.divide(numerator[1], w0_lower)
    target = Decimal(block.k_upper_text)
    passed = k_req_upper < target
    return {
        "y_prime": block.y_prime,
        "x_prime": block.x_prime,
        "target_k_upper_exact_decimal": block.k_upper_text,
        "reciprocal_prime_block_interval": interval_row(reciprocal_block),
        "rho_y_interval": interval_row(rho_y),
        "rho_x_interval": interval_row(rho_x),
        "k_req_numerator_interval": interval_row(numerator),
        "w0_certified_lower": str(w0_lower),
        "k_req_certified_upper": str(k_req_upper),
        "target_minus_k_req_upper": str(target - k_req_upper),
        "status": "PASS" if passed else "FAIL",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Proof-oriented Decimal interval certificate for five finite "
            "CA endpoint margins and four MVDC K_req targets."
        )
    )
    parser.add_argument("--precision", type=int, default=40)
    parser.add_argument("--gamma-terms", type=int, default=100_000)
    parser.add_argument("--series-cutoff", type=int, default=10_000)
    parser.add_argument("--profile-search-limit", type=int, default=100_000)
    parser.add_argument("--w0-steps", type=int, default=5_000)
    parser.add_argument("--progress-every", type=int, default=250_000)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started = time.perf_counter()
    arithmetic = DecimalIntervals(args.precision)
    max_endpoint = max(spec.support_prime for spec in ENDPOINTS)
    sieve_limit = max_endpoint + 2_000

    sieve_started = time.perf_counter()
    sieve = odd_sieve(sieve_limit)
    sieve_seconds = time.perf_counter() - sieve_started
    for spec in ENDPOINTS:
        if not is_prime(sieve, spec.support_prime):
            raise ValueError(f"declared endpoint {spec.support_prime} is not prime")

    search_primes = list(
        primes_from_sieve(sieve, args.profile_search_limit)
    )
    audit_started = time.perf_counter()
    profile_audits = [
        audit_profile(
            arithmetic,
            spec,
            next_prime(sieve, spec.support_prime, sieve_limit),
            search_primes,
        )
        for spec in ENDPOINTS
    ]
    profile_audit_seconds = time.perf_counter() - audit_started

    gamma_started = time.perf_counter()
    gamma = harmonic_gamma_interval(arithmetic, args.gamma_terms)
    gamma_seconds = time.perf_counter() - gamma_started

    endpoint_set = {spec.support_prime for spec in ENDPOINTS}
    theta_checkpoint_set = {
        endpoint
        for audit in profile_audits
        for endpoint in audit.layer_endpoints
    }
    scan_started = time.perf_counter()
    prefix_values, theta_checkpoints, prime_count, prime_digest = (
        scan_prime_prefixes(
            arithmetic,
            sieve,
            max_endpoint,
            endpoint_set,
            theta_checkpoint_set,
            args.series_cutoff,
            args.progress_every,
        )
    )
    prime_scan_seconds = time.perf_counter() - scan_started
    missing = sorted(endpoint_set - prefix_values.keys())
    missing_theta = sorted(theta_checkpoint_set - theta_checkpoints.keys())
    if missing or missing_theta:
        raise ValueError(
            f"missing endpoints={missing}, missing theta checkpoints={missing_theta}"
        )

    endpoint_started = time.perf_counter()
    endpoint_rows: list[dict[str, object]] = []
    gap_by_endpoint: dict[int, Interval] = {}
    for audit in profile_audits:
        row, gap, _ = endpoint_certificate(
            arithmetic,
            audit,
            prefix_values[audit.support_prime],
            theta_checkpoints,
            search_primes,
            gamma,
        )
        endpoint_rows.append(row)
        gap_by_endpoint[audit.support_prime] = gap
    endpoint_seconds = time.perf_counter() - endpoint_started

    block_started = time.perf_counter()
    block_rows = [
        block_certificate(
            arithmetic, block, prefix_values, gap_by_endpoint, args.w0_steps
        )
        for block in BLOCKS
    ]
    block_seconds = time.perf_counter() - block_started
    total_seconds = time.perf_counter() - started

    endpoint_pass = all(row["status"] == "PASS" for row in endpoint_rows)
    profile_pass = all(
        audit.unique_on_open_support_interval for audit in profile_audits
    )
    block_pass = all(row["status"] == "PASS" for row in block_rows)
    all_passed = endpoint_pass and profile_pass and block_pass

    payload: dict[str, object] = {
        "certificate_status": "PASS" if all_passed else "FAIL",
        "method": {
            "endpoint_identity": (
                "G_CA=A_CA+B_log-E equals gamma+log(log(log(n_CA)))"
                "-log(sigma(n_CA)/n_CA). The program encloses this exact log gap."
            ),
            "profile_rule": (
                "Exponent a_p is the number of layers j with epsilon<tau_j(p), "
                "where tau_j(p)=log((1-p^(-j-1))/(1-p^(-j)))/log(p)."
            ),
            "profile_monotonicity": (
                "For fixed j, tau_j(p) decreases with p: its numerator is "
                "log(1+(p-1)/(p*(p^j-1))), a positive decreasing function, "
                "while log(p) increases. Thus two boundary-prime comparisons "
                "certify an entire layer."
            ),
            "gamma_certificate": (
                "The program forms H_n by directed Decimal summation and uses "
                "1/(2n+1)<H_n-log(n)-gamma<1/(2n)."
            ),
            "large_prime_log_sigma": (
                "For p above the configured cutoff, S_2<=log(1+1/p)<=S_3 "
                "from the alternating logarithm series; all terms are rational "
                "and outward rounded."
            ),
            "k_req_identity": (
                "When U=W(Y) and R=W(x)+G(x), K_req's numerator equals "
                "log(log x/log Y)-sum_{Y<p<=x}(1/p)+rho(x)-rho(Y)-G(x). "
                "The second-order product tail cancels exactly."
            ),
            "w0_lower_bound": (
                "After z=log(t), W0 is the integral of "
                "exp(-z/2)*(z+1)/z^2. This integrand is positive and decreasing; "
                "directed right rectangles over an inner endpoint interval give "
                "a certified lower bound."
            ),
        },
        "soundness_assumptions": [
            (
                "CPython Decimal ln, exp, and sqrt obey their documented "
                "correctly-rounded ROUND_HALF_EVEN semantics. Every such result "
                "is padded by one ulp before interval propagation."
            ),
            (
                "Python integer arithmetic and the included odd-only Eratosthenes "
                "sieve execute correctly; no external prime or CA table is trusted."
            ),
            (
                "The elementary harmonic-number inequality and the stated "
                "monotonicity of the CA transition thresholds are used as analytic lemmas."
            ),
        ],
        "scope_and_tie_limitations": [
            "This is a finite certificate for the five displayed supports, not an asymptotic theorem.",
            (
                "The audit proves one exponent profile throughout each open support "
                "parameter interval tau_1(q)<epsilon<tau_1(x), and verifies the "
                "chosen exact decimal epsilon strictly inside it."
            ),
            (
                "At epsilon=tau_1(x), exact support x forces a_x=1; at "
                "epsilon=tau_1(q), exact support x forces a_q=0. The strict "
                "higher-layer separation checks therefore extend profile uniqueness "
                "to the closed support interval for these five supports."
            ),
            (
                "The K_req rows evaluate the exact-reset choice U=W(Y) algebraically; "
                "W(Y) cancels and is not independently enclosed. They certify only "
                "K_req<K_target. A separate signed-cell certificate is needed to "
                "prove K>=K_target. This is not an asymptotic propagation theorem."
            ),
        ],
        "configuration": {
            "precision_decimal_digits": args.precision,
            "gamma_harmonic_terms": args.gamma_terms,
            "large_prime_series_cutoff": args.series_cutoff,
            "profile_search_limit": args.profile_search_limit,
            "w0_right_rectangle_steps": args.w0_steps,
            "max_sieve_endpoint": max_endpoint,
        },
        "runtime_environment": {
            "python_version": sys.version,
            "decimal_version": decimal.__version__,
            "libmpdec_version": decimal.__libmpdec_version__,
            "invocation_argv": sys.argv,
            "certifier_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "euler_gamma_interval": interval_row(gamma),
        "prime_count_through_max_endpoint": prime_count,
        "prime_stream_sha256_little_endian_uint64": prime_digest,
        "endpoints": endpoint_rows,
        "blocks": block_rows,
        "timing_seconds": {
            "sieve": sieve_seconds,
            "profile_audit": profile_audit_seconds,
            "gamma": gamma_seconds,
            "prime_scan": prime_scan_seconds,
            "endpoint_evaluation": endpoint_seconds,
            "block_evaluation": block_seconds,
            "total": total_seconds,
        },
    }
    args.report.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(f"wrote {args.report}")
    print(
        f"status={payload['certificate_status']} primes={prime_count:,} "
        f"scan={prime_scan_seconds:.2f}s total={total_seconds:.2f}s"
    )
    for row in endpoint_rows:
        gap = row["exact_log_robin_gap_interval"]
        print(
            f"x={row['support_prime']:>9} {row['status']} "
            f"G>={gap['lower']} scaled>={row['scaled_gap_interval']['lower']}"
        )
    for row in block_rows:
        print(
            f"{row['y_prime']}->{row['x_prime']} {row['status']} "
            f"K_req<={row['k_req_certified_upper']} "
            f"target={row['target_k_upper_exact_decimal']}"
        )


if __name__ == "__main__":
    main()
