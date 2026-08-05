from __future__ import annotations

import argparse
import bisect
import decimal
import hashlib
import json
import math
import platform
import sys
import time
from dataclasses import dataclass
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_REPORT = ROOT / "ca_exact_buffer_interval_certificate_report.json"

DEFAULT_MIN_PRIME = 3_299
DEFAULT_MAX_PRIME = 56_048_351
DEFAULT_TARGET = "0.7825"
CURRENT_SUPPORTS = {3_329_267, 6_382_007, 12_253_883, 29_093_377, 56_048_351}

Interval = tuple[Decimal, Decimal]


class DecimalIntervals:
    """Outward-rounded intervals over the documented CPython Decimal backend."""

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
        if interval[0] <= 0:
            raise ValueError("logarithm interval must be positive")
        return (
            self._padded_unary("ln", interval[0])[0],
            self._padded_unary("ln", interval[1])[1],
        )

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

    def neg(self, value: Interval) -> Interval:
        return -value[1], -value[0]

    def mul(self, left: Interval, right: Interval) -> Interval:
        lower_values = [self.floor.multiply(a, b) for a in left for b in right]
        upper_values = [self.ceiling.multiply(a, b) for a in left for b in right]
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
        numerator_decimal = Decimal(numerator)
        denominator_decimal = Decimal(denominator)
        return (
            self.floor.divide(numerator_decimal, denominator_decimal),
            self.ceiling.divide(numerator_decimal, denominator_decimal),
        )

    def fixed_interval(
        self, lower_integer: int, upper_integer: int, decimal_digits: int
    ) -> Interval:
        exponent = Decimal(-decimal_digits)
        return (
            self.floor.scaleb(Decimal(lower_integer), exponent),
            self.ceiling.scaleb(Decimal(upper_integer), exponent),
        )


@dataclass(frozen=True)
class OptimizerEvent:
    prime: int
    layer: int
    threshold: Interval
    activation_index: int
    activation_support: int
    delta_a: Interval
    log_prime: Interval
    active_gap: Decimal
    inactive_gap: Decimal | None


@dataclass(frozen=True)
class PrimePowerEvent:
    value: int
    prime: int
    exponent: int


def interval_row(value: Interval) -> dict[str, str]:
    return {"lower": str(value[0]), "upper": str(value[1])}


def odd_sieve(limit: int) -> bytearray:
    """Return a deterministic odd-only Eratosthenes table through limit."""

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


def primes_from_sieve(flags: bytearray, limit: int) -> list[int]:
    return [2] + [
        value for value in range(3, limit + 1, 2) if flags[value // 2]
    ]


def transition_tau(
    arithmetic: DecimalIntervals, prime: int, layer: int
) -> Interval:
    """Enclose the CA threshold for changing the exponent from layer-1 to layer."""

    numerator = prime ** (layer + 1) - 1
    denominator = prime * (prime**layer - 1)
    log_ratio = arithmetic.ln(arithmetic.rational(numerator, denominator))
    return arithmetic.div(log_ratio, arithmetic.ln_integer(prime))


def transition_delta_a(
    arithmetic: DecimalIntervals, prime: int, layer: int
) -> Interval:
    """Enclose the change in A when an exponent rises from layer-1 to layer."""

    numerator = prime ** (layer + 1) - 1
    denominator = prime * (prime**layer - 1)
    return arithmetic.neg(
        arithmetic.ln(arithmetic.rational(numerator, denominator))
    )


def optimizer_control(arithmetic: DecimalIntervals, value: int) -> Interval:
    """Enclose 1/(x log x), the exact separable R_core optimizer threshold."""

    denominator = arithmetic.mul(
        arithmetic.point(value), arithmetic.ln_integer(value)
    )
    return arithmetic.div(arithmetic.point(1), denominator)


def ca_support_tau1(arithmetic: DecimalIntervals, value: int) -> Interval:
    """Enclose tau_1(x), used only to audit the superseded CA-only stream."""

    return arithmetic.div(
        arithmetic.ln(arithmetic.rational(value + 1, value)),
        arithmetic.ln_integer(value),
    )


def prime_power_events(primes: list[int], limit: int) -> list[PrimePowerEvent]:
    events: list[PrimePowerEvent] = []
    for prime in primes:
        if prime * prime > limit:
            break
        exponent = 2
        value = prime * prime
        while value <= limit:
            events.append(PrimePowerEvent(value, prime, exponent))
            if value > limit // prime:
                break
            value *= prime
            exponent += 1
    events.sort(key=lambda event: (event.value, event.prime, event.exponent))
    return events


def locate_activation(
    arithmetic: DecimalIntervals,
    primes: list[int],
    threshold: Interval,
    control_cache: dict[int, Interval],
) -> tuple[int, Decimal, Decimal | None]:
    """Locate the first support with 1/(x log x) strictly below a threshold."""

    def cached_control(index: int) -> Interval:
        value = primes[index]
        if value not in control_cache:
            control_cache[value] = optimizer_control(arithmetic, value)
        return control_cache[value]

    low = 0
    high = len(primes)
    while low < high:
        middle = (low + high) // 2
        support_control = cached_control(middle)
        if support_control[1] < threshold[0]:
            high = middle
        elif support_control[0] > threshold[1]:
            low = middle + 1
        else:
            raise ArithmeticError(
                "unresolved optimizer-threshold comparison at support "
                f"{primes[middle]}"
            )
    if low >= len(primes):
        raise ArithmeticError("event does not activate inside the sieved range")

    active_control = cached_control(low)
    active_gap = arithmetic.floor.subtract(threshold[0], active_control[1])
    if active_gap <= 0:
        raise ArithmeticError("activation-side threshold gap is not positive")

    inactive_gap: Decimal | None = None
    if low:
        inactive_control = cached_control(low - 1)
        inactive_gap = arithmetic.floor.subtract(
            inactive_control[0], threshold[1]
        )
        if inactive_gap <= 0:
            raise ArithmeticError("inactive-side threshold gap is not positive")
    return low, active_gap, inactive_gap


def build_optimizer_events(
    arithmetic: DecimalIntervals,
    primes: list[int],
    max_prime: int,
) -> tuple[list[OptimizerEvent], dict[str, object]]:
    """Build every full-support R_core optimizer event through max_prime."""

    max_control = optimizer_control(arithmetic, max_prime)
    cutoff = math.isqrt(2 * max_prime)
    if (cutoff + 1) * (cutoff + 1) <= 2 * max_prime:
        cutoff += 1
    candidate_end = 0
    while candidate_end < len(primes) and primes[candidate_end] <= cutoff:
        candidate_end += 1

    raw_events: list[tuple[int, int, Interval]] = []
    endpoint_overlaps: list[dict[str, object]] = []
    for prime in primes[:candidate_end]:
        layer = 2
        while True:
            threshold = transition_tau(arithmetic, prime, layer)
            if threshold[0] > max_control[1]:
                raw_events.append((prime, layer, threshold))
                layer += 1
                continue
            if threshold[1] < max_control[0]:
                break
            endpoint_overlaps.append(
                {
                    "prime": prime,
                    "layer": layer,
                    "threshold": interval_row(threshold),
                    "optimizer_control_at_max": interval_row(max_control),
                }
            )
            break

    if endpoint_overlaps:
        raise ArithmeticError(
            f"{len(endpoint_overlaps)} event comparisons overlap 1/(x log x) at max"
        )

    control_cache: dict[int, Interval] = {max_prime: max_control}
    log_cache: dict[int, Interval] = {}
    events: list[OptimizerEvent] = []
    minimum_active_gap: Decimal | None = None
    minimum_inactive_gap: Decimal | None = None
    activation_digest = hashlib.sha256()
    per_prime_last_activation: dict[int, int] = {}

    for prime, layer, threshold in raw_events:
        activation_index, active_gap, inactive_gap = locate_activation(
            arithmetic, primes, threshold, control_cache
        )
        if activation_index < per_prime_last_activation.get(prime, 0):
            raise AssertionError("higher layers activated out of order")
        per_prime_last_activation[prime] = activation_index
        if primes[activation_index] < prime:
            raise AssertionError("an optimizer event activates before its own prime")
        if prime not in log_cache:
            log_cache[prime] = arithmetic.ln_integer(prime)
        event = OptimizerEvent(
            prime=prime,
            layer=layer,
            threshold=threshold,
            activation_index=activation_index,
            activation_support=primes[activation_index],
            delta_a=transition_delta_a(arithmetic, prime, layer),
            log_prime=log_cache[prime],
            active_gap=active_gap,
            inactive_gap=inactive_gap,
        )
        events.append(event)
        minimum_active_gap = (
            active_gap
            if minimum_active_gap is None
            else min(minimum_active_gap, active_gap)
        )
        if inactive_gap is not None:
            minimum_inactive_gap = (
                inactive_gap
                if minimum_inactive_gap is None
                else min(minimum_inactive_gap, inactive_gap)
            )
        activation_digest.update(prime.to_bytes(8, "little"))
        activation_digest.update(layer.to_bytes(4, "little"))
        activation_digest.update(primes[activation_index].to_bytes(8, "little"))

    events.sort(key=lambda event: (event.activation_index, event.prime, event.layer))
    return events, {
        "candidate_prime_cutoff": cutoff,
        "candidate_primes_checked": candidate_end,
        "event_count_through_max": len(events),
        "threshold_endpoint_overlap_count": len(endpoint_overlaps),
        "minimum_certified_active_side_gap": str(minimum_active_gap),
        "minimum_certified_inactive_side_gap": str(minimum_inactive_gap),
        "activation_stream_sha256": activation_digest.hexdigest(),
        "control_cache_size_during_event_location": len(control_cache),
    }


def exact_small_prime_entry(
    arithmetic: DecimalIntervals, prime: int
) -> Interval:
    """Enclose 1/p-log(1+1/p), the net insertion into A-T."""

    return arithmetic.sub(
        arithmetic.rational(1, prime),
        arithmetic.ln(arithmetic.rational(prime + 1, prime)),
    )


def fixed_series_entry(
    prime: int, scale: int
) -> tuple[int, int]:
    """Bound 1/p-log(1+1/p) by the S5/S6 alternating partial sums."""

    p2 = prime * prime
    p3 = p2 * prime
    p4 = p3 * prime
    p5 = p4 * prime
    p6 = p5 * prime
    lower_numerator = 30 * p3 - 20 * p2 + 15 * prime - 12
    lower_denominator = 60 * p5
    upper_numerator = 30 * p4 - 20 * p3 + 15 * p2 - 12 * prime + 10
    upper_denominator = 60 * p6
    lower_integer = scale * lower_numerator // lower_denominator
    upper_integer = (
        scale * upper_numerator + upper_denominator - 1
    ) // upper_denominator
    return lower_integer, upper_integer


def scaled_buffer(
    arithmetic: DecimalIntervals,
    support: int,
    a_minus_t: Interval,
    h_minus_l: Interval,
) -> tuple[Interval, Interval, Interval]:
    sqrt_x = arithmetic.sqrt_integer(support)
    log_x = arithmetic.ln_integer(support)
    scale = arithmetic.mul(sqrt_x, log_x)
    first = arithmetic.mul(a_minus_t, scale)
    second = arithmetic.div(h_minus_l, sqrt_x)
    return arithmetic.add(first, second), sqrt_x, log_x


def direct_state_check(
    arithmetic: DecimalIntervals,
    primes: list[int],
    support_index: int,
    optimizer_events: list[OptimizerEvent],
) -> dict[str, object]:
    """Independently reconstruct A-T and H-L at one small checkpoint."""

    support = primes[support_index]
    exponent_by_prime = {prime: 1 for prime in primes[: support_index + 1]}
    for event in optimizer_events:
        if event.activation_index > support_index:
            break
        exponent_by_prime[event.prime] += 1

    a_core = arithmetic.point(0)
    h_value = arithmetic.point(0)
    l_value = arithmetic.point(0)
    t_value = arithmetic.point(0)
    for prime in primes[: support_index + 1]:
        exponent = exponent_by_prime[prime]
        factor = arithmetic.rational(
            prime ** (exponent + 1) - 1,
            prime ** (exponent + 1),
        )
        a_core = arithmetic.sub(a_core, arithmetic.ln(factor))
        if exponent > 1:
            h_value = arithmetic.add(
                h_value,
                arithmetic.mul(
                    arithmetic.point(exponent - 1),
                    arithmetic.ln_integer(prime),
                ),
            )

        higher_tail = arithmetic.sub(
            arithmetic.neg(
                arithmetic.ln(arithmetic.rational(prime - 1, prime))
            ),
            arithmetic.rational(1, prime),
        )
        t_value = arithmetic.add(t_value, higher_tail)
        value = prime * prime
        power = 2
        while value <= support:
            l_value = arithmetic.add(l_value, arithmetic.ln_integer(prime))
            t_value = arithmetic.sub(
                t_value, arithmetic.rational(1, power * value)
            )
            value *= prime
            power += 1

    direct_a_minus_t = arithmetic.sub(a_core, t_value)
    direct_h_minus_l = arithmetic.sub(h_value, l_value)
    direct_scaled, _sqrt_x, _log_x = scaled_buffer(
        arithmetic, support, direct_a_minus_t, direct_h_minus_l
    )
    return {
        "support_prime": support,
        "direct_a_minus_t": interval_row(direct_a_minus_t),
        "direct_h_minus_l": interval_row(direct_h_minus_l),
        "direct_scaled_buffer": interval_row(direct_scaled),
        "profile_max_exponent": max(exponent_by_prime.values()),
    }


def scan(args: argparse.Namespace) -> dict[str, object]:
    arithmetic = DecimalIntervals(args.precision)
    started = time.perf_counter()
    sieve_started = time.perf_counter()
    sieve = odd_sieve(args.max_prime)
    primes = primes_from_sieve(sieve, args.max_prime)
    sieve_seconds = time.perf_counter() - sieve_started
    if args.min_prime not in primes:
        raise ValueError("min-prime must itself be prime for exact range coverage")
    if args.max_prime not in primes:
        raise ValueError("max-prime must itself be prime for exact range coverage")

    first_index = primes.index(args.min_prime)
    last_index = len(primes) - 1
    if primes[last_index] != args.max_prime:
        raise AssertionError("prime stream does not end at max-prime")

    event_started = time.perf_counter()
    optimizer_events, optimizer_metadata = build_optimizer_events(
        arithmetic, primes, args.max_prime
    )
    power_events = prime_power_events(primes, args.max_prime)
    event_seconds = time.perf_counter() - event_started

    optimizer_buckets: dict[int, list[OptimizerEvent]] = {}
    for event in optimizer_events:
        optimizer_buckets.setdefault(event.activation_index, []).append(event)

    power_index = 0
    small_entry = arithmetic.point(0)
    other_a_minus_t = arithmetic.point(0)
    h_minus_l = arithmetic.point(0)
    fixed_lower = 0
    fixed_upper = 0
    fixed_scale = 10**args.fixed_scale_digits
    target = Decimal(args.target)
    prime_digest = hashlib.sha256()
    power_digest = hashlib.sha256()
    minimum_lower: Decimal | None = None
    minimum_row: dict[str, object] | None = None
    maximum_interval_width = Decimal(0)
    pass_count = 0
    fail_count = 0
    inconclusive_count = 0
    rows: list[dict[str, object]] = []
    active_optimizer_events = 0
    max_bucket_size = 0
    scan_started = time.perf_counter()

    requested_direct_supports = [
        int(part.strip())
        for part in args.direct_check_primes.split(",")
        if part.strip()
    ]
    if not requested_direct_supports:
        raise ValueError("direct-check-primes must contain at least one integer")
    direct_supports = set()
    for requested in requested_direct_supports:
        direct_support_index = bisect.bisect_right(primes, requested) - 1
        if direct_support_index < 0:
            raise ValueError("a direct-check-prime is below the first sieved prime")
        direct_supports.add(primes[direct_support_index])
    checkpoint_set = set(CURRENT_SUPPORTS)
    checkpoint_set.update({args.min_prime, args.max_prime})
    checkpoint_set.update(direct_supports)

    for index, prime in enumerate(primes):
        prime_digest.update(prime.to_bytes(8, "little"))
        if prime <= args.series_cutoff:
            small_entry = arithmetic.add(
                small_entry, exact_small_prime_entry(arithmetic, prime)
            )
        else:
            lower_integer, upper_integer = fixed_series_entry(prime, fixed_scale)
            fixed_lower += lower_integer
            fixed_upper += upper_integer

        while power_index < len(power_events) and power_events[power_index].value <= prime:
            event = power_events[power_index]
            other_a_minus_t = arithmetic.add(
                other_a_minus_t,
                arithmetic.rational(1, event.exponent * event.value),
            )
            h_minus_l = arithmetic.sub(
                h_minus_l, arithmetic.ln_integer(event.prime)
            )
            power_digest.update(event.value.to_bytes(8, "little"))
            power_digest.update(event.prime.to_bytes(8, "little"))
            power_digest.update(event.exponent.to_bytes(4, "little"))
            power_index += 1

        bucket = optimizer_buckets.get(index, [])
        max_bucket_size = max(max_bucket_size, len(bucket))
        for event in bucket:
            other_a_minus_t = arithmetic.add(other_a_minus_t, event.delta_a)
            h_minus_l = arithmetic.add(h_minus_l, event.log_prime)
            active_optimizer_events += 1

        if index < first_index:
            continue

        fixed_entry = arithmetic.fixed_interval(
            fixed_lower, fixed_upper, args.fixed_scale_digits
        )
        a_minus_t = arithmetic.add(
            arithmetic.add(small_entry, fixed_entry), other_a_minus_t
        )
        value, sqrt_x, log_x = scaled_buffer(
            arithmetic, prime, a_minus_t, h_minus_l
        )
        width = arithmetic.ceiling.subtract(value[1], value[0])
        lower_margin = arithmetic.floor.subtract(value[0], target)
        maximum_interval_width = max(maximum_interval_width, width)
        if value[0] > target:
            pass_count += 1
            classification = "PASS"
        elif value[1] <= target:
            fail_count += 1
            classification = "FAIL"
        else:
            inconclusive_count += 1
            classification = "INCONCLUSIVE"

        if minimum_lower is None or value[0] < minimum_lower:
            minimum_lower = value[0]
            minimum_row = {
                "support_prime": prime,
                "scaled_buffer": interval_row(value),
                "lower_margin_over_target": str(lower_margin),
                "a_minus_t": interval_row(a_minus_t),
                "h_minus_l": interval_row(h_minus_l),
                "sqrt_x": interval_row(sqrt_x),
                "log_x": interval_row(log_x),
                "active_optimizer_events": active_optimizer_events,
                "processed_prime_power_events": power_index,
                "classification": classification,
            }

        if prime in checkpoint_set:
            rows.append(
                {
                    "support_prime": prime,
                    "scaled_buffer": interval_row(value),
                    "lower_margin_over_target": str(lower_margin),
                    "a_minus_t": interval_row(a_minus_t),
                    "h_minus_l": interval_row(h_minus_l),
                    "active_optimizer_events": active_optimizer_events,
                    "processed_prime_power_events": power_index,
                    "classification": classification,
                }
            )

        processed = index - first_index + 1
        if args.progress_every and processed % args.progress_every == 0:
            elapsed = time.perf_counter() - scan_started
            print(
                f"certified {processed:,} supports through {prime:,} "
                f"in {elapsed:.1f}s; minimum lower={minimum_lower}",
                flush=True,
            )

    scan_seconds = time.perf_counter() - scan_started
    if power_index != len(power_events):
        raise AssertionError("not all prime-power events were processed")
    if active_optimizer_events != len(optimizer_events):
        raise AssertionError("not all optimizer events were processed")
    assert minimum_row is not None and minimum_lower is not None

    direct_checks = []
    checked_supports = sorted(
        {args.min_prime}
        | {support for support in direct_supports if support <= args.max_prime}
    )
    for support in checked_supports:
        check_index = bisect.bisect_right(primes, support) - 1
        check = direct_state_check(arithmetic, primes, check_index, optimizer_events)
        stream_row = next(
            (row for row in rows if row["support_prime"] == primes[check_index]), None
        )
        if stream_row is None:
            raise AssertionError("direct reconstruction checkpoint is missing")
        overlap_checks: dict[str, bool] = {}
        for direct_key, stream_key in (
            ("direct_a_minus_t", "a_minus_t"),
            ("direct_h_minus_l", "h_minus_l"),
            ("direct_scaled_buffer", "scaled_buffer"),
        ):
            direct_interval = check[direct_key]
            stream_interval = stream_row[stream_key]
            overlap_checks[direct_key] = not (
                Decimal(direct_interval["upper"])
                < Decimal(stream_interval["lower"])
                or Decimal(stream_interval["upper"])
                < Decimal(direct_interval["lower"])
            )
        if not all(overlap_checks.values()):
            raise AssertionError(
                f"direct reconstruction does not overlap stream at {support}"
            )
        check["stream_checkpoint"] = stream_row
        check["interval_overlap_checks"] = overlap_checks
        direct_checks.append(check)

    fixed_interval = arithmetic.fixed_interval(
        fixed_lower, fixed_upper, args.fixed_scale_digits
    )
    status = (
        "PASS"
        if fail_count == 0 and inconclusive_count == 0
        else "FAIL"
        if fail_count > 0
        else "INCONCLUSIVE"
    )
    return {
        "claim": {
            "quantity": (
                "sqrt(x) log(x) D_*(x), where D_*(x)=R_*(x)-C_pp(x)"
            ),
            "range": (
                f"every prime support {args.min_prime} <= x <= {args.max_prime}"
            ),
            "strict_target": args.target,
            "optimizer": (
                "R_*(x) is the minimum of R_core=A+H/(x log x) over all "
                "full-support exponent vectors a_p>=1 for p<=x"
            ),
        },
        "status": status,
        "status_counts": {
            "support_count": last_index - first_index + 1,
            "pass": pass_count,
            "fail": fail_count,
            "inconclusive": inconclusive_count,
        },
        "minimum": minimum_row,
        "maximum_interval_width": str(maximum_interval_width),
        "checkpoints": sorted(rows, key=lambda row: int(row["support_prime"])),
        "optimizer_event_certificate": {
            **optimizer_metadata,
            "events_active_at_first_support": sum(
                event.activation_index <= first_index for event in optimizer_events
            ),
            "events_activating_inside_target_range": sum(
                first_index < event.activation_index <= last_index
                for event in optimizer_events
            ),
            "maximum_simultaneous_activation_bucket": max_bucket_size,
            "tie_rule": (
                "At tau_s(p)=1/(x log x), the adjacent exponents give exactly "
                "the same R_core contribution; excluding the event is tie-safe."
            ),
            "minimum_argument": (
                "R_core is separable over primes. Raising an exponent from s-1 "
                "to s changes it by log(p)*(1/(x log x)-tau_s(p)). Therefore "
                "the transition is used exactly when tau_s(p)>1/(x log x); at "
                "equality either adjacent exponent is minimizing."
            ),
            "exact_optimizer_formula": (
                "a_p=1+#{s>=2: tau_s(p)>1/(x log x)} for p<=x; equality may "
                "be included or excluded"
            ),
        },
        "correction_audit": {
            "superseded_condition": "tau_s(p)>tau_1(x)",
            "superseded_scope": (
                "This gives the minimum only inside the family of CA profiles "
                "having exact support x. It is not the full-support R_* optimizer."
            ),
            "correct_condition": "tau_s(p)>1/(x log x)",
            "effect": (
                "The report and PASS status use only the corrected full-support "
                "optimizer stream. No result from the superseded full run is used."
            ),
        },
        "prime_power_certificate": {
            "event_count": len(power_events),
            "event_stream_sha256": power_digest.hexdigest(),
            "identity": (
                "C_pp=(psi-theta)/(x log x)+"
                "sum_{p<=x,m>=2,p^m>x}1/(m p^m)"
            ),
        },
        "combined_stream_identity": {
            "state": ["A_minus_T=A-T", "H_minus_L=H-(psi-theta)"],
            "scaled_buffer": (
                "(A-T)sqrt(x)log(x)+(H-(psi-theta))/sqrt(x)"
            ),
            "prime_insertion": (
                "1/p-log(1+1/p), bounded by exact logarithms for small p "
                "and alternating S5/S6 series for large p"
            ),
            "optimizer_transition": (
                "A changes by log(1-p^-s)-log(1-p^(-s-1)); H changes by log p"
            ),
            "prime_power_transition": (
                "A-T changes by 1/(m p^m); H-(psi-theta) changes by -log p"
            ),
        },
        "series_certificate": {
            "exact_log_cutoff": args.series_cutoff,
            "fixed_scale_decimal_digits": args.fixed_scale_digits,
            "large_prime_lower_partial_sum": (
                "1/(2p^2)-1/(3p^3)+1/(4p^4)-1/(5p^5)"
            ),
            "large_prime_upper_partial_sum": "lower+1/(6p^6)",
            "final_large_prime_sum_interval": interval_row(fixed_interval),
        },
        "direct_reconstruction_checks": direct_checks,
        "coverage": {
            "first_support": primes[first_index],
            "last_support": primes[last_index],
            "prime_count_through_max": len(primes),
            "prime_stream_sha256_little_endian_uint64": prime_digest.hexdigest(),
        },
        "rigor_scope": {
            "analytic_lemmas_used": [
                "1/(x log x) is strictly decreasing for x>1",
                "tau_s(p) is strictly decreasing in p and in s",
                (
                    "events with p>sqrt(2*max_prime) are inactive: "
                    "p^2>2X and log(p)>log(X)/2 imply p^2 log(p)>X log(X), "
                    "hence tau_s(p)<=tau_2(p)<1/(p^2 log(p))<1/(X log(X)); "
                    "the control is still larger for every x<X"
                ),
                (
                    "the alternating logarithm terms decrease for 0<1/p<1, "
                    "so S5 is a lower bound and S6 is an upper bound"
                ),
            ],
            "computational_assumptions": [
                (
                    "CPython Decimal ln and sqrt obey their documented "
                    "correct-rounding semantics; one ulp is added on both sides"
                ),
                "directed Decimal arithmetic has the requested IEEE rounding mode",
                "the deterministic Eratosthenes sieve implementation is correct",
                "Python arbitrary-precision integer arithmetic is exact",
            ],
            "formalization_gap": (
                "This is an executable interval certificate, not a proof-assistant "
                "formalization. The four elementary analytic lemmas and the runtime "
                "semantics above remain trusted inputs."
            ),
        },
        "parameters": {
            "min_prime": args.min_prime,
            "max_prime": args.max_prime,
            "target": args.target,
            "precision": args.precision,
            "series_cutoff": args.series_cutoff,
            "fixed_scale_digits": args.fixed_scale_digits,
            "direct_check_primes": requested_direct_supports,
        },
        "timing_seconds": {
            "sieve": sieve_seconds,
            "event_construction": event_seconds,
            "support_scan": scan_seconds,
            "total_before_report_serialization": time.perf_counter() - started,
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Certify the exact full-support finite CA core buffer with "
            "outward-rounded intervals."
        )
    )
    parser.add_argument("--min-prime", type=int, default=DEFAULT_MIN_PRIME)
    parser.add_argument("--max-prime", type=int, default=DEFAULT_MAX_PRIME)
    parser.add_argument("--target", default=DEFAULT_TARGET)
    parser.add_argument("--precision", type=int, default=40)
    parser.add_argument("--series-cutoff", type=int, default=1_000)
    parser.add_argument("--fixed-scale-digits", type=int, default=70)
    parser.add_argument(
        "--direct-check-primes",
        default="5297,5303,10007",
        help=(
            "Comma-separated direct reconstruction checkpoints. The default "
            "straddles a known CA-only versus full-support optimizer shift."
        ),
    )
    parser.add_argument("--progress-every", type=int, default=250_000)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.min_prime < 3:
        raise ValueError("min-prime must be at least 3")
    if args.max_prime < args.min_prime:
        raise ValueError("max-prime must be at least min-prime")
    if args.series_cutoff < 2:
        raise ValueError("series-cutoff must be at least 2")
    if args.fixed_scale_digits < args.precision + 10:
        raise ValueError("fixed-scale-digits must exceed precision by at least 10")

    report = scan(args)
    script_path = Path(__file__).resolve()
    report["runtime_environment"] = {
        "python": sys.version,
        "platform": platform.platform(),
        "decimal_module": decimal.__file__,
        "libmpdec_version": getattr(decimal, "__libmpdec_version__", "unknown"),
        "argv": sys.argv,
        "certifier_script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    minimum = report["minimum"]
    print(f"wrote {args.report}")
    print(
        f"status={report['status']} supports={report['status_counts']['support_count']:,} "
        f"minimum_lower={minimum['scaled_buffer']['lower']} "
        f"at x={minimum['support_prime']} "
        f"margin={minimum['lower_margin_over_target']}"
    )


if __name__ == "__main__":
    main()
