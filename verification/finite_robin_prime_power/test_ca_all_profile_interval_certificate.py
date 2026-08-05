from __future__ import annotations

import argparse
import bisect
import decimal
import hashlib
import json
import math
import sys
import time
from array import array
from dataclasses import dataclass
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN
from pathlib import Path
from typing import Iterator


ROOT = Path(__file__).resolve().parent
DEFAULT_REPORT = ROOT / "ca_all_profile_interval_certificate_report.json"
DEFAULT_ENDPOINT_CROSSCHECK = ROOT / "ca_endpoint_interval_certificate_report.json"
CURRENT_SUPPORT_CHECKPOINTS = (
    3_329_267,
    6_382_007,
    12_253_883,
    29_093_377,
    56_048_351,
)

Interval = tuple[Decimal, Decimal]


class DecimalIntervals:
    """Minimal outward-rounded Decimal interval arithmetic."""

    def __init__(self, precision: int) -> None:
        if precision < 30:
            raise ValueError("precision must be at least 30 decimal digits")
        self.precision = precision
        self.nearest = Context(prec=precision, rounding=ROUND_HALF_EVEN)
        self.floor = Context(prec=precision, rounding=ROUND_FLOOR)
        self.ceiling = Context(prec=precision, rounding=ROUND_CEILING)
        self.padding = Context(prec=precision + 8, rounding=ROUND_HALF_EVEN)
        self._ulp_cache: dict[int, Decimal] = {}
        self.log1p_even_order_counts = {2: 0, 4: 0}

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

    def ln(self, value: Interval) -> Interval:
        if value[0] <= 0:
            raise ValueError("logarithm interval must be positive")
        return (
            self._padded_unary("ln", value[0])[0],
            self._padded_unary("ln", value[1])[1],
        )

    def ln_lower(self, value: Decimal) -> Decimal:
        if value <= 0:
            raise ValueError("logarithm input must be positive")
        return self._padded_unary("ln", value)[0]

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
        if left[0] >= 0 and right[0] >= 0:
            return (
                self.floor.multiply(left[0], right[0]),
                self.ceiling.multiply(left[1], right[1]),
            )
        if left[1] <= 0 and right[1] <= 0:
            return (
                self.floor.multiply(left[1], right[1]),
                self.ceiling.multiply(left[0], right[0]),
            )
        if left[0] >= 0 and right[1] <= 0:
            return (
                self.floor.multiply(left[1], right[0]),
                self.ceiling.multiply(left[0], right[1]),
            )
        if left[1] <= 0 and right[0] >= 0:
            return (
                self.floor.multiply(left[0], right[1]),
                self.ceiling.multiply(left[1], right[0]),
            )
        lower_values = [self.floor.multiply(a, b) for a in left for b in right]
        upper_values = [self.ceiling.multiply(a, b) for a in left for b in right]
        return min(lower_values), max(upper_values)

    def div(self, numerator: Interval, denominator: Interval) -> Interval:
        if denominator[0] <= 0 <= denominator[1]:
            raise ZeroDivisionError("division interval contains zero")
        if denominator[0] > 0 and numerator[0] >= 0:
            return (
                self.floor.divide(numerator[0], denominator[1]),
                self.ceiling.divide(numerator[1], denominator[0]),
            )
        if denominator[0] > 0 and numerator[1] <= 0:
            return (
                self.floor.divide(numerator[0], denominator[0]),
                self.ceiling.divide(numerator[1], denominator[1]),
            )
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


@dataclass(frozen=True)
class TransitionEvent:
    prime: int
    layer: int
    tau: Interval
    log_prime: Interval
    log_sigma_increment: Interval
    approximate_tau: float


@dataclass
class ProfileState:
    log_n: Interval
    log_sigma_over_n: Interval
    log_log_n: Interval | None = None
    gap: Interval | None = None


def interval_row(value: Interval) -> dict[str, str]:
    return {"lower": str(value[0]), "upper": str(value[1])}


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


def primes_from_sieve(flags: bytearray, limit: int) -> Iterator[int]:
    yield 2
    for value in range(3, limit + 1, 2):
        if flags[value // 2]:
            yield value


def build_prime_array(flags: bytearray, limit: int) -> array:
    result = array("I")
    result.extend(primes_from_sieve(flags, limit))
    return result


def transition_data(
    arithmetic: DecimalIntervals, prime: int, layer: int
) -> tuple[Interval, Interval, Interval]:
    """Return log(p), the local log-sigma gain, and its CA threshold."""

    numerator = prime ** (layer + 1) - 1
    denominator = prime * (prime**layer - 1)
    log_prime = arithmetic.ln_integer(prime)
    log_gain = arithmetic.ln(arithmetic.rational(numerator, denominator))
    tau = arithmetic.div(log_gain, log_prime)
    return log_prime, log_gain, tau


def transition_tau_float(prime: int, layer: int) -> float:
    inverse_power = prime ** (-layer)
    log_gain = math.log1p(-inverse_power / prime) - math.log1p(-inverse_power)
    return log_gain / math.log(prime)


def squarefree_sigma_log_factor(
    arithmetic: DecimalIntervals, prime: int, series_cutoff: int
) -> Interval:
    """Enclose log(1+1/p), with an alternating rational tail."""

    if prime <= series_cutoff:
        return arithmetic.ln(arithmetic.rational(prime + 1, prime))
    first = arithmetic.rational(1, prime)
    second = arithmetic.rational(1, 2 * prime * prime)
    third = arithmetic.rational(1, 3 * prime * prime * prime)
    lower = arithmetic.sub(first, second)[0]
    upper = arithmetic.add(arithmetic.sub(first, second), third)[1]
    return lower, upper


def harmonic_gamma_interval(
    arithmetic: DecimalIntervals, term_count: int
) -> Interval:
    """Self-certify Euler's constant using harmonic-number inequalities."""

    harmonic = arithmetic.point(0)
    for denominator in range(1, term_count + 1):
        harmonic = arithmetic.add(harmonic, arithmetic.rational(1, denominator))
    log_n = arithmetic.ln_integer(term_count)
    lower = arithmetic.sub(
        arithmetic.sub(harmonic, log_n),
        arithmetic.rational(1, 2 * term_count),
    )[0]
    upper = arithmetic.sub(
        arithmetic.sub(harmonic, log_n),
        arithmetic.rational(1, 2 * term_count + 1),
    )[1]
    return lower, upper


def gap_lower_bound(
    arithmetic: DecimalIntervals, state: ProfileState, gamma: Interval
) -> Decimal:
    """One-sided enclosure of gamma+log(log(log(n)))-log(sigma(n)/n)."""

    if state.gap is not None:
        return state.gap[0]
    log_log_n_lower = arithmetic.ln_lower(state.log_n[0])
    log_log_log_n_lower = arithmetic.ln_lower(log_log_n_lower)
    partial = arithmetic.floor.add(gamma[0], log_log_log_n_lower)
    return arithmetic.floor.subtract(partial, state.log_sigma_over_n[1])


def full_gap_interval(
    arithmetic: DecimalIntervals, state: ProfileState, gamma: Interval
) -> Interval:
    log_log_n = arithmetic.ln(state.log_n)
    log_log_log_n = arithmetic.ln(log_log_n)
    return arithmetic.sub(arithmetic.add(gamma, log_log_log_n), state.log_sigma_over_n)


def log1p_small_positive_interval(
    arithmetic: DecimalIntervals, value: Interval
) -> Interval:
    """Enclose log(1+z) between consecutive even/odd partial sums."""

    if value[0] < 0 or value[1] >= 1:
        raise ValueError("the small positive log1p enclosure requires 0<=z<1")
    even_order = 2 if value[1] <= Decimal("0.00001") else 4
    arithmetic.log1p_even_order_counts[even_order] += 1
    odd_order = even_order + 1
    power = value
    partial = arithmetic.point(0)
    even_partial: Interval | None = None
    odd_partial: Interval | None = None
    for exponent in range(1, odd_order + 1):
        term = arithmetic.div(power, arithmetic.point(exponent))
        partial = (
            arithmetic.add(partial, term)
            if exponent % 2 == 1
            else arithmetic.sub(partial, term)
        )
        if exponent == even_order:
            even_partial = partial
        elif exponent == odd_order:
            odd_partial = partial
        if exponent < odd_order:
            power = arithmetic.mul(power, value)
    assert even_partial is not None and odd_partial is not None
    return even_partial[0], odd_partial[1]


def add_profile_increment(
    arithmetic: DecimalIntervals,
    state: ProfileState,
    log_n_increment: Interval,
    log_sigma_increment: Interval,
) -> None:
    """Update the profile and, after initialization, its gap without new ln calls."""

    if state.log_log_n is not None or state.gap is not None:
        if state.log_log_n is None or state.gap is None:
            raise ValueError("derived profile intervals must be initialized together")
        relative_log_n_increment = arithmetic.div(
            log_n_increment, state.log_n
        )
        log_log_n_increment = log1p_small_positive_interval(
            arithmetic, relative_log_n_increment
        )
        relative_log_log_increment = arithmetic.div(
            log_log_n_increment, state.log_log_n
        )
        log_log_log_increment = log1p_small_positive_interval(
            arithmetic, relative_log_log_increment
        )
        state.log_log_n = arithmetic.add(
            state.log_log_n, log_log_n_increment
        )
        state.gap = arithmetic.add(
            state.gap,
            arithmetic.sub(log_log_log_increment, log_sigma_increment),
        )
    state.log_n = arithmetic.add(state.log_n, log_n_increment)
    state.log_sigma_over_n = arithmetic.add(
        state.log_sigma_over_n, log_sigma_increment
    )


def add_transition(
    arithmetic: DecimalIntervals, state: ProfileState, event: TransitionEvent
) -> None:
    add_profile_increment(
        arithmetic, state, event.log_prime, event.log_sigma_increment
    )


def add_first_layer(
    arithmetic: DecimalIntervals,
    state: ProfileState,
    prime: int,
    series_cutoff: int,
) -> Interval:
    log_prime = arithmetic.ln_integer(prime)
    add_profile_increment(
        arithmetic,
        state,
        log_prime,
        squarefree_sigma_log_factor(arithmetic, prime, series_cutoff),
    )
    return log_prime


def find_support_bucket_float(
    primes: array,
    min_index: int,
    max_index: int,
    event_tau: float,
) -> int:
    """Locate i with tau_1(p_i)>event_tau>tau_1(p_{i+1})."""

    low = min_index
    high = max_index + 2
    while low < high:
        middle = (low + high) // 2
        if transition_tau_float(primes[middle], 1) > event_tau:
            low = middle + 1
        else:
            high = middle
    return low - 1


def generate_higher_events(
    arithmetic: DecimalIntervals,
    primes: array,
    min_index: int,
    max_index: int,
    higher_search_limit: int,
) -> tuple[
    list[TransitionEvent],
    list[TransitionEvent],
    list[dict[str, object]],
    dict[str, object],
]:
    """Generate every layer >=2 event relevant to the support sweep."""

    upper_boundary = transition_data(arithmetic, primes[min_index], 1)[2]
    lower_boundary = transition_data(arithmetic, primes[max_index + 1], 1)[2]
    initial_events: list[TransitionEvent] = []
    active_events: list[TransitionEvent] = []
    unresolved: list[dict[str, object]] = []
    largest_layer = 1
    primes_examined = 0
    certified_layer_stops = 0

    search_end = bisect.bisect_right(primes, higher_search_limit)
    for prime in primes[:search_end]:
        primes_examined += 1
        layer = 2
        while True:
            log_prime, log_gain, tau = transition_data(arithmetic, prime, layer)
            largest_layer = max(largest_layer, layer)
            if tau[1] < lower_boundary[0]:
                certified_layer_stops += 1
                break
            event = TransitionEvent(
                prime=prime,
                layer=layer,
                tau=tau,
                log_prime=log_prime,
                log_sigma_increment=log_gain,
                approximate_tau=transition_tau_float(prime, layer),
            )
            if tau[0] > upper_boundary[1]:
                initial_events.append(event)
            elif tau[1] < upper_boundary[0]:
                active_events.append(event)
            else:
                unresolved.append(
                    {
                        "type": "global_upper_boundary_overlap",
                        "prime": prime,
                        "layer": layer,
                        "event_tau": interval_row(tau),
                        "boundary_tau": interval_row(upper_boundary),
                    }
                )
            layer += 1

    # For s>=2, tau_s(p)<=tau_2(p)<1/(p(p+1)log(p)).  The value at the
    # configured integer limit therefore excludes every prime beyond it.
    log_limit = arithmetic.ln_integer(higher_search_limit)
    tail_bound = arithmetic.div(
        arithmetic.point(1),
        arithmetic.mul(
            arithmetic.point(higher_search_limit * (higher_search_limit + 1)),
            log_limit,
        ),
    )
    tail_excluded = tail_bound[1] < lower_boundary[0]
    coverage = {
        "global_upper_tau1_min_support": interval_row(upper_boundary),
        "global_lower_tau1_next_after_max_support": interval_row(lower_boundary),
        "higher_prime_search_limit": higher_search_limit,
        "primes_examined_for_higher_layers": primes_examined,
        "largest_layer_examined_including_first_below_range": largest_layer,
        "certified_layer_stop_count": certified_layer_stops,
        "tail_tau_upper_bound": interval_row(tail_bound),
        "all_higher_events_beyond_search_limit_excluded": tail_excluded,
    }
    return initial_events, active_events, unresolved, coverage


def assign_event_buckets(
    arithmetic: DecimalIntervals,
    primes: array,
    min_index: int,
    max_index: int,
    events: list[TransitionEvent],
) -> tuple[dict[int, list[TransitionEvent]], list[dict[str, object]]]:
    buckets: dict[int, list[TransitionEvent]] = {}
    unresolved: list[dict[str, object]] = []

    for event in events:
        candidate = find_support_bucket_float(
            primes, min_index, max_index, event.approximate_tau
        )
        certified_index: int | None = None
        for index in range(max(min_index, candidate - 3), min(max_index, candidate + 3) + 1):
            upper = transition_data(arithmetic, primes[index], 1)[2]
            lower = transition_data(arithmetic, primes[index + 1], 1)[2]
            if upper[0] > event.tau[1] and event.tau[0] > lower[1]:
                certified_index = index
                break
        if certified_index is None:
            nearest = max(min_index, min(max_index, candidate))
            unresolved.append(
                {
                    "type": "support_boundary_or_assignment_overlap",
                    "event_prime": event.prime,
                    "event_layer": event.layer,
                    "event_tau": interval_row(event.tau),
                    "candidate_support_prime": primes[nearest],
                    "candidate_upper_tau": interval_row(
                        transition_data(arithmetic, primes[nearest], 1)[2]
                    ),
                    "candidate_lower_tau": interval_row(
                        transition_data(arithmetic, primes[nearest + 1], 1)[2]
                    ),
                }
            )
            continue
        buckets.setdefault(certified_index, []).append(event)

    for index, bucket in buckets.items():
        bucket.sort(key=lambda item: item.approximate_tau, reverse=True)
        for earlier, later in zip(bucket, bucket[1:]):
            if earlier.tau[0] <= later.tau[1]:
                unresolved.append(
                    {
                        "type": "simultaneous_or_unresolved_higher_tie_group",
                        "support_prime": primes[index],
                        "first_event": {
                            "prime": earlier.prime,
                            "layer": earlier.layer,
                            "tau": interval_row(earlier.tau),
                        },
                        "second_event": {
                            "prime": later.prime,
                            "layer": later.layer,
                            "tau": interval_row(later.tau),
                        },
                    }
                )
    return buckets, unresolved


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Certify the exact logarithmic Robin gap for every CA exponent "
            "profile, including both sides of every higher-layer tie."
        )
    )
    parser.add_argument("--min-support", type=int, default=3_329_267)
    parser.add_argument("--max-support", type=int, default=56_048_351)
    parser.add_argument("--precision", type=int, default=44)
    parser.add_argument("--gamma-terms", type=int, default=200_000)
    parser.add_argument("--series-cutoff", type=int, default=20_000)
    parser.add_argument("--higher-search-limit", type=int, default=100_000)
    parser.add_argument("--summary-block-size", type=int, default=250_000)
    parser.add_argument("--progress-every", type=int, default=250_000)
    parser.add_argument("--sieve-padding", type=int, default=1_000)
    parser.add_argument(
        "--endpoint-crosscheck-report",
        type=Path,
        default=DEFAULT_ENDPOINT_CROSSCHECK,
    )
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    started = time.perf_counter()
    arithmetic = DecimalIntervals(args.precision)
    if args.min_support > args.max_support:
        raise ValueError("min-support must not exceed max-support")
    if args.summary_block_size <= 0:
        raise ValueError("summary-block-size must be positive")

    sieve_limit = args.max_support + args.sieve_padding
    flags = odd_sieve(sieve_limit)
    primes = build_prime_array(flags, sieve_limit)
    min_index = bisect.bisect_left(primes, args.min_support)
    max_index = bisect.bisect_left(primes, args.max_support)
    if min_index >= len(primes) or primes[min_index] != args.min_support:
        raise ValueError(f"min-support {args.min_support} is not prime")
    if max_index >= len(primes) or primes[max_index] != args.max_support:
        raise ValueError(f"max-support {args.max_support} is not prime")
    if max_index + 1 >= len(primes):
        raise ValueError("sieve padding did not reach the next prime")
    if args.higher_search_limit >= len(flags) * 2:
        raise ValueError("higher-search-limit lies outside the sieve")

    gamma = harmonic_gamma_interval(arithmetic, args.gamma_terms)
    initial_events, active_events, unresolved, higher_coverage = generate_higher_events(
        arithmetic,
        primes,
        min_index,
        max_index,
        args.higher_search_limit,
    )
    buckets, bucket_unresolved = assign_event_buckets(
        arithmetic, primes, min_index, max_index, active_events
    )
    unresolved.extend(bucket_unresolved)

    event_digest = hashlib.sha256()
    event_structure_digest = hashlib.sha256()
    for event in sorted(initial_events, key=lambda item: (item.prime, item.layer)):
        event_digest.update(b"I")
        event_structure_digest.update(b"I")
        event_digest.update(event.prime.to_bytes(8, byteorder="little", signed=False))
        event_structure_digest.update(
            event.prime.to_bytes(8, byteorder="little", signed=False)
        )
        event_digest.update(event.layer.to_bytes(4, byteorder="little", signed=False))
        event_structure_digest.update(
            event.layer.to_bytes(4, byteorder="little", signed=False)
        )
        event_digest.update(f"{event.tau[0]}|{event.tau[1]}".encode("ascii"))
    for support_index in sorted(buckets):
        for event in buckets[support_index]:
            event_digest.update(b"A")
            event_structure_digest.update(b"A")
            event_digest.update(
                primes[support_index].to_bytes(8, byteorder="little", signed=False)
            )
            event_structure_digest.update(
                primes[support_index].to_bytes(
                    8, byteorder="little", signed=False
                )
            )
            event_digest.update(
                event.prime.to_bytes(8, byteorder="little", signed=False)
            )
            event_structure_digest.update(
                event.prime.to_bytes(8, byteorder="little", signed=False)
            )
            event_digest.update(
                event.layer.to_bytes(4, byteorder="little", signed=False)
            )
            event_structure_digest.update(
                event.layer.to_bytes(4, byteorder="little", signed=False)
            )
            event_digest.update(f"{event.tau[0]}|{event.tau[1]}".encode("ascii"))
    higher_event_stream_sha256 = event_digest.hexdigest()
    higher_event_structure_sha256 = event_structure_digest.hexdigest()

    classification_pass = (
        not unresolved
        and bool(higher_coverage["all_higher_events_beyond_search_limit_excluded"])
    )
    if not classification_pass:
        report = {
            "certificate_status": "FAIL",
            "failure_stage": "higher_transition_classification",
            "configuration": {
                key: str(value) if isinstance(value, Path) else value
                for key, value in vars(args).items()
            },
            "higher_transition_coverage": higher_coverage,
            "higher_event_stream_sha256": higher_event_stream_sha256,
            "higher_event_structure_sha256": higher_event_structure_sha256,
            "unresolved_transition_relations": unresolved,
            "timing_seconds": str(Decimal(str(time.perf_counter() - started))),
        }
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(json.dumps(report, indent=2))
        raise SystemExit(1)

    state = ProfileState(arithmetic.point(0), arithmetic.point(0))
    prime_digest = hashlib.sha256()
    log_prime_at_min: Interval | None = None
    prefix_started = time.perf_counter()
    for index in range(min_index + 1):
        prime = primes[index]
        log_prime_at_min = add_first_layer(
            arithmetic, state, prime, args.series_cutoff
        )
        prime_digest.update(prime.to_bytes(8, byteorder="little", signed=False))
        if args.progress_every and (index + 1) % args.progress_every == 0:
            print(
                f"prefix: {index + 1:,} primes through {prime:,} "
                f"in {time.perf_counter() - prefix_started:.1f}s",
                flush=True,
            )
    assert log_prime_at_min is not None
    for event in initial_events:
        add_transition(arithmetic, state, event)
    state.log_log_n = arithmetic.ln(state.log_n)
    state.gap = full_gap_interval(arithmetic, state, gamma)

    support_count = max_index - min_index + 1
    evaluated_count = 0
    passed_count = 0
    failed_count = 0
    failure_examples: list[dict[str, object]] = []
    transition_ties: list[dict[str, object]] = []
    block_summaries: list[dict[str, object]] = []
    checkpoint_snapshots: dict[int, list[dict[str, object]]] = {
        support: []
        for support in CURRENT_SUPPORT_CHECKPOINTS
        if args.min_support <= support <= args.max_support
    }
    global_minimum: dict[str, object] | None = None
    block_minimum: dict[str, object] | None = None
    maximum_interval_widths: dict[str, dict[str, object]] = {}
    block_first_index = min_index
    last_gap_lower: Decimal | None = None

    def evaluate(
        support_index: int,
        profile_kind: str,
        triggering_event: TransitionEvent | None,
    ) -> Decimal:
        nonlocal evaluated_count, passed_count, failed_count
        nonlocal global_minimum, block_minimum
        support = primes[support_index]
        lower = gap_lower_bound(arithmetic, state, gamma)
        width_inputs = {
            "log_n": state.log_n,
            "log_sigma_over_n": state.log_sigma_over_n,
            "log_log_n": state.log_log_n,
            "gap": state.gap,
        }
        for quantity, value in width_inputs.items():
            assert value is not None
            width = arithmetic.ceiling.subtract(value[1], value[0])
            previous = maximum_interval_widths.get(quantity)
            if previous is None or width > Decimal(str(previous["width"])):
                maximum_interval_widths[quantity] = {
                    "width": str(width),
                    "support_prime": support,
                    "profile_kind": profile_kind,
                    "triggering_higher_event": (
                        None
                        if triggering_event is None
                        else {
                            "prime": triggering_event.prime,
                            "layer": triggering_event.layer,
                        }
                    ),
                }
        evaluated_count += 1
        if lower > 0:
            passed_count += 1
        else:
            failed_count += 1
            if len(failure_examples) < 100:
                failure_examples.append(
                    {
                        "support_prime": support,
                        "profile_kind": profile_kind,
                        "event_prime": None if triggering_event is None else triggering_event.prime,
                        "event_layer": None if triggering_event is None else triggering_event.layer,
                        "gap_lower": str(lower),
                    }
                )
        candidate = {
            "support_prime": support,
            "support_prime_index_zero_based": support_index,
            "profile_kind": profile_kind,
            "triggering_higher_event": (
                None
                if triggering_event is None
                else {"prime": triggering_event.prime, "layer": triggering_event.layer}
            ),
            "gap_lower": str(lower),
            "state": ProfileState(
                state.log_n,
                state.log_sigma_over_n,
                state.log_log_n,
                state.gap,
            ),
        }
        if support in checkpoint_snapshots:
            checkpoint_snapshots[support].append(
                {
                    "profile_kind": profile_kind,
                    "triggering_higher_event": candidate[
                        "triggering_higher_event"
                    ],
                    "gap_lower": str(lower),
                    "state": ProfileState(
                        state.log_n,
                        state.log_sigma_over_n,
                        state.log_log_n,
                        state.gap,
                    ),
                }
            )
        if global_minimum is None or lower < Decimal(str(global_minimum["gap_lower"])):
            global_minimum = candidate
        if block_minimum is None or lower < Decimal(str(block_minimum["gap_lower"])):
            block_minimum = candidate
        return lower

    current_log_prime = log_prime_at_min
    scan_started = time.perf_counter()
    last_gap_lower = evaluate(min_index, "after_forced_first_layer", None)

    for support_index in range(min_index, max_index + 1):
        support = primes[support_index]
        if support_index != min_index:
            # The first-layer state was evaluated immediately after the preceding
            # boundary; its value remains the pre-event state here.
            assert last_gap_lower is not None

        for event in buckets.get(support_index, []):
            pre_lower = last_gap_lower
            add_transition(arithmetic, state, event)
            post_lower = evaluate(
                support_index, "after_higher_transition", event
            )
            transition_ties.append(
                {
                    "support_prime": support,
                    "event_prime": event.prime,
                    "event_layer": event.layer,
                    "event_tau_interval": interval_row(event.tau),
                    "tie_branch_exponent_before": event.layer - 1,
                    "tie_branch_exponent_after": event.layer,
                    "gap_lower_before_branch": str(pre_lower),
                    "gap_lower_after_branch": str(post_lower),
                    "both_tie_branches_certified_positive": (
                        pre_lower > 0 and post_lower > 0
                    ),
                }
            )
            last_gap_lower = post_lower

        is_block_end = (
            (support_index - min_index + 1) % args.summary_block_size == 0
            or support_index == max_index
        )
        if is_block_end:
            assert block_minimum is not None
            block_summaries.append(
                {
                    "first_support_prime": primes[block_first_index],
                    "last_support_prime": support,
                    "support_count": support_index - block_first_index + 1,
                    "minimum_profile_gap_lower": block_minimum["gap_lower"],
                    "minimum_at_support_prime": block_minimum["support_prime"],
                    "minimum_profile_kind": block_minimum["profile_kind"],
                    "minimum_triggering_higher_event": block_minimum[
                        "triggering_higher_event"
                    ],
                }
            )
            block_first_index = support_index + 1
            block_minimum = None

        if support_index < max_index:
            next_support = primes[support_index + 1]
            current_log_prime = add_first_layer(
                arithmetic, state, next_support, args.series_cutoff
            )
            prime_digest.update(
                next_support.to_bytes(8, byteorder="little", signed=False)
            )
            last_gap_lower = evaluate(
                support_index + 1, "after_forced_first_layer", None
            )

        processed_supports = support_index - min_index + 1
        if args.progress_every and processed_supports % args.progress_every == 0:
            print(
                f"profiles: {processed_supports:,}/{support_count:,} supports; "
                f"{evaluated_count:,} profiles; min lower "
                f"{global_minimum['gap_lower'] if global_minimum else 'n/a'}; "
                f"{time.perf_counter() - scan_started:.1f}s",
                flush=True,
            )

    assert global_minimum is not None
    minimum_state = global_minimum.pop("state")
    assert isinstance(minimum_state, ProfileState)
    assert minimum_state.gap is not None
    minimum_interval = minimum_state.gap
    minimum_direct_interval = full_gap_interval(
        arithmetic, minimum_state, gamma
    )
    minimum_support = int(global_minimum["support_prime"])
    minimum_log_prime = arithmetic.ln_integer(minimum_support)
    scaled_minimum = arithmetic.mul(
        minimum_interval,
        arithmetic.mul(
            arithmetic.sqrt_integer(minimum_support), minimum_log_prime
        ),
    )
    scaled_minimum_direct = arithmetic.mul(
        minimum_direct_interval,
        arithmetic.mul(
            arithmetic.sqrt_integer(minimum_support), minimum_log_prime
        ),
    )
    all_ties_pass = all(
        bool(row["both_tie_branches_certified_positive"])
        for row in transition_ties
    )
    profile_count_identity = evaluated_count == support_count + len(active_events)
    maximum_gap_width = Decimal(str(maximum_interval_widths["gap"]["width"]))
    interval_width_control_pass = maximum_gap_width < Decimal(
        str(global_minimum["gap_lower"])
    )

    endpoint_crosscheck_path = args.endpoint_crosscheck_report.resolve()
    endpoint_crosscheck_available = endpoint_crosscheck_path.is_file()
    old_endpoint_by_support: dict[int, dict[str, object]] = {}
    endpoint_crosscheck_sha256: str | None = None
    if endpoint_crosscheck_available:
        endpoint_crosscheck_bytes = endpoint_crosscheck_path.read_bytes()
        endpoint_crosscheck_sha256 = hashlib.sha256(
            endpoint_crosscheck_bytes
        ).hexdigest()
        endpoint_crosscheck_data = json.loads(endpoint_crosscheck_bytes)
        old_endpoint_by_support = {
            int(row["support_prime"]): row
            for row in endpoint_crosscheck_data.get("endpoints", [])
        }

    checkpoint_rows: list[dict[str, object]] = []
    checkpoint_crosschecks_pass = endpoint_crosscheck_available
    for support in sorted(checkpoint_snapshots):
        support_index = bisect.bisect_left(primes, support)
        profile_rows: list[dict[str, object]] = []
        for snapshot in checkpoint_snapshots[support]:
            snapshot_state = snapshot.pop("state")
            assert isinstance(snapshot_state, ProfileState)
            assert snapshot_state.gap is not None
            streamed_gap_interval = snapshot_state.gap
            direct_gap_interval = full_gap_interval(
                arithmetic, snapshot_state, gamma
            )
            scale = arithmetic.mul(
                arithmetic.sqrt_integer(support),
                arithmetic.ln_integer(support),
            )
            profile_rows.append(
                {
                    **snapshot,
                    "streamed_gap_interval": interval_row(
                        streamed_gap_interval
                    ),
                    "direct_decimal_gap_interval": interval_row(
                        direct_gap_interval
                    ),
                    "streamed_scaled_gap_interval_sqrt_x_log_x": interval_row(
                        arithmetic.mul(streamed_gap_interval, scale)
                    ),
                    "direct_scaled_gap_interval_sqrt_x_log_x": interval_row(
                        arithmetic.mul(direct_gap_interval, scale)
                    ),
                }
            )

        old_row = old_endpoint_by_support.get(support)
        old_interval: Interval | None = None
        if old_row is not None:
            old_gap = old_row["exact_log_robin_gap_interval"]
            assert isinstance(old_gap, dict)
            old_interval = (Decimal(old_gap["lower"]), Decimal(old_gap["upper"]))
        unique_profile = len(profile_rows) == 1
        streamed_overlap: bool | None = None
        direct_overlap: bool | None = None
        if unique_profile and old_interval is not None:
            streamed_gap = profile_rows[0]["streamed_gap_interval"]
            direct_gap = profile_rows[0]["direct_decimal_gap_interval"]
            assert isinstance(streamed_gap, dict) and isinstance(direct_gap, dict)
            streamed_overlap = (
                Decimal(streamed_gap["lower"]) <= old_interval[1]
                and old_interval[0] <= Decimal(streamed_gap["upper"])
            )
            direct_overlap = (
                Decimal(direct_gap["lower"]) <= old_interval[1]
                and old_interval[0] <= Decimal(direct_gap["upper"])
            )
        row_pass = (
            unique_profile
            and streamed_overlap is True
            and direct_overlap is True
        )
        checkpoint_crosschecks_pass = checkpoint_crosschecks_pass and row_pass
        checkpoint_rows.append(
            {
                "support_prime": support,
                "higher_transition_event_count_inside_support_interval": len(
                    buckets.get(support_index, [])
                ),
                "numerical_profile_count": len(profile_rows),
                "expected_unique_profile_for_this_current_support": True,
                "profiles": profile_rows,
                "existing_endpoint_certificate_gap_interval": (
                    None if old_interval is None else interval_row(old_interval)
                ),
                "streamed_interval_overlap_with_existing_endpoint_certificate": streamed_overlap,
                "direct_interval_overlap_with_existing_endpoint_certificate": direct_overlap,
                "checkpoint_crosscheck_status": "PASS" if row_pass else "FAIL",
            }
        )

    overall_pass = (
        classification_pass
        and failed_count == 0
        and passed_count == evaluated_count
        and all_ties_pass
        and profile_count_identity
        and checkpoint_crosschecks_pass
        and interval_width_control_pass
    )

    script_path = Path(__file__).resolve()
    report = {
        "certificate_status": "PASS" if overall_pass else "FAIL",
        "claim_certified": (
            "For every prime support x in the configured inclusive range, "
            "the exact logarithmic Robin gap G(n)=gamma+log(log(log(n)))"
            "-log(sigma(n)/n) is positive for every CA exponent profile "
            "compatible with exact support x, including both exponent choices "
            "at every higher-layer transition tie."
        ),
        "method": {
            "ca_profile_partition": (
                "The CA objective is additive over prime exponents.  Adding the "
                "s-th factor of p has gain log((1-p^(-s-1))/(1-p^(-s)))"
                "-epsilon*log(p), so its threshold is tau_s(p).  The sweep adds "
                "first layers at consecutive prime-support boundaries and all "
                "strictly interlaced higher-layer events in decreasing threshold "
                "order.  Each event's pre- and post-state are the two CA profiles "
                "at the exact tie."
            ),
            "tie_completeness": (
                "Interval comparisons certify that every active higher event lies "
                "strictly inside one support interval, that adjacent active events "
                "are distinct, and hence every tie group has size one in this run. "
                "At a first-layer boundary, exact support forces the old next prime "
                "to be absent and the new support prime to be present; both forced "
                "states are evaluated by the adjacent sweep rows."
            ),
            "higher_event_exhaustion": (
                "For fixed p, tau_s decreases with s because the exact gain ratio "
                "is 1+(p-1)/(p*(p^s-1)).  For s>=2 and p beyond the search limit, "
                "tau_s(p)<=tau_2(p)<1/(p*(p+1)*log(p)); the reported outward "
                "interval upper bound is below the global lower epsilon boundary."
            ),
            "gap_arithmetic": (
                "log(n) and log(sigma(n)/n) are maintained incrementally with "
                "separate outward-rounded Decimal intervals.  After a direct "
                "Decimal initialization, log(log(n)) and G are updated exactly "
                "through two positive log1p increments.  Each log1p is enclosed "
                "by the explicit alternating-series parity bracket S_(2m)"
                "<=log1p<=S_(2m+1), using m=1 for z<=1e-5 and m=2 "
                "otherwise.  The lower endpoint of the streamed G interval is "
                "checked at every numerical profile."
            ),
        },
        "soundness_notes": [
            "CPython Decimal ln and sqrt are correctly rounded under ROUND_HALF_EVEN; one ulp is subtracted/added after every transcendental call.",
            "All rational operations and accumulated sums use directed ROUND_FLOOR/ROUND_CEILING contexts.",
            "Euler's constant is enclosed using 1/(2N+1)<H_N-log(N)-gamma<1/(2N).",
            "For primes above the series cutoff, z-z^2/2<=log(1+z)<=z-z^2/2+z^3/3 with z=1/p.",
            "Every streamed log1p argument is checked to lie in [0,1), and an explicit consecutive even/odd partial-sum remainder bracket is applied.",
        ],
        "configuration": {
            "min_support_prime": args.min_support,
            "max_support_prime": args.max_support,
            "next_prime_after_max_support": primes[max_index + 1],
            "precision_decimal_digits": args.precision,
            "gamma_harmonic_terms": args.gamma_terms,
            "large_prime_series_cutoff": args.series_cutoff,
            "higher_prime_search_limit": args.higher_search_limit,
            "summary_block_size": args.summary_block_size,
            "sieve_limit": sieve_limit,
            "endpoint_crosscheck_report": str(endpoint_crosscheck_path),
        },
        "runtime_environment": {
            "python_version": sys.version,
            "decimal_module_version": getattr(decimal, "__version__", None),
            "libmpdec_version": getattr(decimal, "__libmpdec_version__", None),
            "argv": sys.argv,
            "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
        },
        "euler_gamma_interval": interval_row(gamma),
        "higher_transition_coverage": higher_coverage,
        "classification": {
            "status": "PASS" if classification_pass else "FAIL",
            "initial_higher_events_strictly_above_range": len(initial_events),
            "active_higher_transition_ties": len(active_events),
            "unresolved_transition_relations": unresolved,
            "simultaneous_higher_tie_groups": 0,
            "every_active_tie_has_two_positive_branches": all_ties_pass,
            "higher_event_stream_sha256": higher_event_stream_sha256,
            "higher_event_structure_sha256": higher_event_structure_sha256,
        },
        "profile_coverage": {
            "prime_support_count": support_count,
            "evaluated_distinct_numerical_profile_count": evaluated_count,
            "expected_count_supports_plus_higher_events": support_count
            + len(active_events),
            "count_identity_pass": profile_count_identity,
            "positive_profile_count": passed_count,
            "nonpositive_or_unresolved_profile_count": failed_count,
            "higher_transition_tie_branch_incidence_count": 2
            * len(active_events),
            "first_layer_internal_boundary_count": support_count - 1,
        },
        "global_uniform_lower_bound": {
            **global_minimum,
            "streamed_gap_interval_at_minimum_lower_bound_profile": interval_row(
                minimum_interval
            ),
            "direct_decimal_gap_interval_at_same_profile": interval_row(
                minimum_direct_interval
            ),
            "streamed_scaled_gap_interval_sqrt_x_log_x": interval_row(
                scaled_minimum
            ),
            "direct_scaled_gap_interval_sqrt_x_log_x": interval_row(
                scaled_minimum_direct
            ),
        },
        "interval_width_audit": {
            "status": "PASS" if interval_width_control_pass else "FAIL",
            "criterion": (
                "The maximum streamed G interval width over every evaluated "
                "support/event profile is strictly smaller than the certified "
                "global positive lower bound."
            ),
            "maximum_widths_over_all_evaluated_profiles": maximum_interval_widths,
            "maximum_gap_width": str(maximum_gap_width),
            "global_positive_lower_bound": global_minimum["gap_lower"],
            "streamed_log1p_even_order_call_counts": {
                str(order): count
                for order, count in arithmetic.log1p_even_order_counts.items()
            },
        },
        "current_support_checkpoints": {
            "status": "PASS" if checkpoint_crosschecks_pass else "FAIL",
            "range_starts_at_exact_prime": 3_329_267,
            "existing_endpoint_certificate_path": str(
                endpoint_crosscheck_path
            ),
            "existing_endpoint_certificate_available": endpoint_crosscheck_available,
            "existing_endpoint_certificate_sha256": endpoint_crosscheck_sha256,
            "rows": checkpoint_rows,
        },
        "support_block_summaries": block_summaries,
        "higher_transition_ties": transition_ties,
        "failure_examples": failure_examples,
        "prime_stream_sha256_little_endian_uint64_through_max_support": prime_digest.hexdigest(),
        "timing_seconds": str(Decimal(str(time.perf_counter() - started))),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(
        f"{report['certificate_status']}: {support_count:,} supports, "
        f"{evaluated_count:,} profiles, {len(active_events):,} higher ties; "
        f"uniform lower bound {global_minimum['gap_lower']}; "
        f"report {args.report}"
    )
    if not overall_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
