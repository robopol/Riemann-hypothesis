from __future__ import annotations

import argparse
import decimal
import hashlib
import json
import platform
import sys
from pathlib import Path

from test_ca_all_profile_interval_certificate import (
    DecimalIntervals,
    interval_row,
    transition_data,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_REPORT = ROOT / "ca_support_11_overlap_certificate_report.json"
INTERVAL_BACKEND = ROOT / "test_ca_all_profile_interval_certificate.py"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Certify that the exact-support-11 CA bridge begins at 55440, "
            "the endpoint of the direct integer base interval."
        )
    )
    parser.add_argument("--precision", type=int, default=60)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    arithmetic = DecimalIntervals(args.precision)
    boundary = transition_data(arithmetic, 11, 1)[2]
    next_boundary = transition_data(arithmetic, 13, 1)[2]

    exponents: dict[int, int] = {}
    comparisons: list[dict[str, object]] = []
    unresolved: list[dict[str, object]] = []
    for prime in (2, 3, 5, 7):
        exponent = 1
        layer = 2
        while True:
            tau = transition_data(arithmetic, prime, layer)[2]
            if tau[0] > boundary[1]:
                relation = "strictly_above_support_11_boundary"
                exponent = layer
            elif tau[1] < boundary[0]:
                relation = "strictly_below_support_11_boundary"
            else:
                relation = "unresolved_overlap"
                unresolved.append(
                    {
                        "prime": prime,
                        "layer": layer,
                        "tau_interval": interval_row(tau),
                    }
                )
            comparisons.append(
                {
                    "prime": prime,
                    "layer": layer,
                    "tau_interval": interval_row(tau),
                    "relation": relation,
                }
            )
            if relation != "strictly_above_support_11_boundary":
                break
            layer += 1
        exponents[prime] = exponent

    # Exact support 11 forces the new first layer to be present at its tie.
    exponents[11] = 1
    ca_value = 1
    previous_value = 1
    for prime, exponent in exponents.items():
        ca_value *= prime**exponent
        if prime != 11:
            previous_value *= prime**exponent

    expected_exponents = {2: 4, 3: 2, 5: 1, 7: 1, 11: 1}
    next_boundary_strictly_lower = next_boundary[1] < boundary[0]
    overall_pass = (
        not unresolved
        and exponents == expected_exponents
        and previous_value == 5_040
        and ca_value == 55_440
        and next_boundary_strictly_lower
    )

    script_path = Path(__file__).resolve()
    report = {
        "certificate_status": "PASS" if overall_pass else "FAIL",
        "claim_certified": (
            "At the first-layer CA threshold tau_1(11), the old profile is "
            "5040 and exact support 11 forces the new profile "
            "2^4*3^2*5*7*11=55440. Thus the support-11 CA bridge overlaps "
            "the direct integer base at 55440."
        ),
        "support_11_boundary_tau_interval": interval_row(boundary),
        "support_13_boundary_tau_interval": interval_row(next_boundary),
        "support_13_boundary_is_strictly_lower": next_boundary_strictly_lower,
        "higher_layer_comparisons": comparisons,
        "unresolved_comparisons": unresolved,
        "old_tie_profile": {
            "integer": previous_value,
            "factorization": "2^4 3^2 5^1 7^1",
            "largest_prime_factor": 7,
        },
        "forced_exact_support_11_profile": {
            "integer": ca_value,
            "factorization": "2^4 3^2 5^1 7^1 11^1",
            "largest_prime_factor": 11,
            "exponents": {str(prime): exponent for prime, exponent in exponents.items()},
        },
        "method": [
            "The CA objective is additive over prime-exponent layers; layer s of p is active above its threshold tau_s(p).",
            "Every displayed higher-layer threshold is interval-separated from tau_1(11).",
            "At the first-layer tie, the old profile excludes 11 and the exact-support-11 profile includes it. No numerical ordering assumption is used at the tie.",
            "The first-layer threshold tau_1(p) decreases with p, and tau_1(13)<tau_1(11), so primes above 11 are absent on the support-11 side.",
        ],
        "runtime_environment": {
            "python_version": sys.version,
            "python_implementation": platform.python_implementation(),
            "decimal_module_version": decimal.__version__,
            "libmpdec_version": decimal.__libmpdec_version__,
            "precision_decimal_digits": args.precision,
            "script_path": str(script_path),
            "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
            "interval_backend_path": str(INTERVAL_BACKEND.resolve()),
            "interval_backend_sha256": hashlib.sha256(
                INTERVAL_BACKEND.read_bytes()
            ).hexdigest(),
            "argv": sys.argv,
        },
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(
        f"{report['certificate_status']}: old={previous_value}, "
        f"exact-support-11={ca_value}; report {args.report.resolve()}"
    )
    if not overall_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
