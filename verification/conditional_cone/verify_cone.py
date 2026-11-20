"""Certify positivity and the strict cone at all prime-power events up to 60,000."""

import argparse
import json
from pathlib import Path
import sys

# The companion archive already contains this standard-library interval backend.
BACKEND = Path(__file__).resolve().parents[1] / "finite_robin_prime_power"
sys.path.insert(0, str(BACKEND))
from test_ca_residual_finite_certificate import (
    DecimalIntervals, euler_gamma_interval, machin_pi_interval,
)


ENDPOINT = 60_000
GAMMA_TERMS = 100_000
PRECISION = 50
REPORT = Path(__file__).with_name("report.json")


def prime_power_bases(limit):
    """Return the base p for each q=p^a, and zero at other indices."""
    composite = bytearray(limit + 1)
    bases = [0] * (limit + 1)
    for p in range(2, limit + 1):
        if composite[p]:
            continue
        for multiple in range(p * p, limit + 1, p):
            composite[multiple] = 1
        power = p
        while power <= limit:
            bases[power] = p
            power *= p
    return bases


def enclose_state(a, q, psi, reciprocal, gamma, log_two_pi, origin):
    """Use inclusive sums for the right derivative immediately after q."""
    x = a.point(q)
    root = a.sqrt(x)
    half = a.rational(1, 2)
    atanh_inverse = a.mul(half, a.ln(a.rational(q + 1, q - 1)))
    u = a.add(a.sub(reciprocal, a.ln(x)),
              a.add(gamma, a.sub(a.rational(1, q), atanh_inverse)))
    v = a.add(a.sub(psi, x),
              a.add(log_two_pi, a.mul(half, a.ln(a.rational(q*q - 1, q*q)))))
    root_u = a.mul(root, u)
    v_over_root = a.div(v, root)
    phi = a.add(a.sub(origin, root_u), v_over_root)
    w = a.mul(half, a.add(root_u, v_over_root))
    k = a.sub(a.sub(root, a.div(a.point(1), a.mul(root, a.point(q*q - 1)))),
              a.mul(a.rational(1, 4), origin))
    radius_squared = a.mul(phi, a.add(a.mul(a.point(2), k), a.mul(a.rational(1, 4), phi)))
    return phi, a.sub(radius_squared, a.mul(w, w))


def verify(precision=PRECISION):
    a = DecimalIntervals(precision)
    gamma = euler_gamma_interval(a, GAMMA_TERMS)
    pi = machin_pi_interval(a)
    log_two_pi = a.ln(a.mul(a.point(2), pi))
    origin = a.sub(a.add(a.point(2), gamma), a.ln(a.mul(a.point(4), pi)))
    psi = reciprocal = a.point(0)
    event_count = 0
    minima = {}
    for q, base in enumerate(prime_power_bases(ENDPOINT)):
        if base == 0:
            continue
        event_count += 1
        weight = a.ln(a.point(base))
        psi = a.add(psi, weight)
        reciprocal = a.add(reciprocal, a.div(weight, a.point(q)))
        phi, determinant = enclose_state(a, q, psi, reciprocal, gamma, log_two_pi, origin)
        for name, interval in (("phi", phi), ("determinant", determinant)):
            lower = interval[0]
            if lower <= 0:
                raise ArithmeticError(f"Failed to certify {name} at q={q}: {interval}")
            if name not in minima or lower < minima[name][0]:
                minima[name] = lower, q
    if event_count != 6148 or a.exp(a.point(11))[1] >= ENDPOINT:
        raise ArithmeticError("The finite range does not match the manuscript")
    return {
        "certified": True, "integer_endpoint": ENDPOINT, "post_event_count": event_count,
        "decimal_precision": precision, "gamma_harmonic_terms": GAMMA_TERMS,
        "minima_of_post_event_lower_endpoints": {
            name: {"lower": str(value), "q": q} for name, (value, q) in minima.items()
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--check", action="store_true", help="compare with the archived report")
    action.add_argument("--write-report", action="store_true", help="regenerate the archived report")
    options = parser.parse_args()
    result = verify()
    if options.check:
        if result != json.loads(REPORT.read_text(encoding="utf-8")):
            raise SystemExit("The recomputed certificate differs from the archived report")
    rendered = json.dumps(result, indent=2) + "\n"
    if options.write_report:
        REPORT.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
