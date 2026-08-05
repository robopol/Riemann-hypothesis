from __future__ import annotations

import argparse
import decimal
import hashlib
import json
import platform
import re
import sys
from dataclasses import dataclass
from decimal import Decimal, ROUND_FLOOR
from pathlib import Path
from typing import Any

from test_ca_all_profile_interval_certificate import DecimalIntervals, interval_row


ROOT = Path(__file__).resolve().parent
DEFAULT_REPORT = ROOT / "ca_all_integer_cutoff_certificate_report.json"
PROFILE_INTERVAL_BACKEND = ROOT / "test_ca_all_profile_interval_certificate.py"
RESIDUAL_INTERVAL_BACKEND = ROOT / "test_ca_residual_finite_certificate.py"
PREVIOUS_UPWARD_CERTIFIER = ROOT / "test_ca_residual_upward_1e21_certificate.py"
DSTAR_CERTIFIER = ROOT / "test_asymptotic_dstar_lower_bound.py"
DSTAR_REPORT = ROOT / "asymptotic_dstar_lower_bound_report.json"

DUSART_THETA_CONSTANT = Decimal("3.965")
SUPPORT_FRACTION_NUMERATOR = 999
SUPPORT_FRACTION_DENOMINATOR = 1_000
CLEAN_DECIMAL_EXPONENT_DIGITS = 2


@dataclass(frozen=True)
class Dependency:
    name: str
    report: Path
    script: Path


DEPENDENCIES = (
    Dependency(
        "integer_base_5041_55440",
        ROOT / "robin_base_5041_55440_interval_report.json",
        ROOT / "test_robin_base_5041_55440_interval.py",
    ),
    Dependency(
        "ca_support_11_overlap_at_55440",
        ROOT / "ca_support_11_overlap_certificate_report.json",
        ROOT / "test_ca_support_11_overlap_certificate.py",
    ),
    Dependency(
        "ca_bridge_support_11_3299",
        ROOT / "ca_all_profile_bridge_11_3299_report.json",
        ROOT / "test_ca_all_profile_interval_certificate.py",
    ),
    Dependency(
        "ca_direct_support_3299_56048351",
        ROOT / "ca_all_profile_interval_certificate_report.json",
        ROOT / "test_ca_all_profile_interval_certificate.py",
    ),
    Dependency(
        "ca_analytic_support_56048351_to_X",
        ROOT / "ca_residual_upward_extended_certificate_report.json",
        ROOT / "test_ca_residual_upward_extended_certificate.py",
    ),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_dependency(
    dependency: Dependency, expected_x_max: int
) -> tuple[dict[str, Any], dict[str, Any]]:
    report = json.loads(dependency.report.read_bytes())
    script_hash = sha256(dependency.script)
    checks: list[dict[str, object]] = []
    transitive_files: list[dict[str, str]] = []

    def record(name: str, passed: bool, observed: object, expected: object) -> None:
        checks.append(
            {
                "name": name,
                "status": "PASS" if passed else "FAIL",
                "observed": observed,
                "expected": expected,
            }
        )

    def record_transitive_hash(
        name: str, embedded_hash: object, path: Path
    ) -> None:
        actual_hash = sha256(path)
        record(name, embedded_hash == actual_hash, embedded_hash, actual_hash)
        transitive_files.append(
            {
                "name": name,
                "path": str(path.resolve()),
                "sha256": actual_hash,
            }
        )

    if dependency.name == "integer_base_5041_55440":
        configuration = report.get("configuration", {})
        coverage = report.get("coverage", {})
        runtime = report.get("runtime_environment", {})
        record("certificate status", report.get("certificate_status") == "PASS", report.get("certificate_status"), "PASS")
        record("base start", configuration.get("start") == 5_041, configuration.get("start"), 5_041)
        record("base end", configuration.get("end") == 55_440, configuration.get("end"), 55_440)
        record("integer count", coverage.get("integer_count") == 50_400, coverage.get("integer_count"), 50_400)
        record("failure count", coverage.get("failure_count") == 0, coverage.get("failure_count"), 0)
        embedded_hash = runtime.get("script_sha256")
        record_transitive_hash(
            "embedded interval backend hash",
            runtime.get("interval_backend_sha256"),
            PROFILE_INTERVAL_BACKEND,
        )
    elif dependency.name == "ca_support_11_overlap_at_55440":
        runtime = report.get("runtime_environment", {})
        old_profile = report.get("old_tie_profile", {})
        new_profile = report.get("forced_exact_support_11_profile", {})
        record("certificate status", report.get("certificate_status") == "PASS", report.get("certificate_status"), "PASS")
        record("old tie profile", old_profile.get("integer") == 5_040, old_profile.get("integer"), 5_040)
        record("forced support-11 profile", new_profile.get("integer") == 55_440, new_profile.get("integer"), 55_440)
        record("forced largest prime", new_profile.get("largest_prime_factor") == 11, new_profile.get("largest_prime_factor"), 11)
        record("unresolved comparisons", report.get("unresolved_comparisons") == [], report.get("unresolved_comparisons"), [])
        embedded_hash = runtime.get("script_sha256")
        record_transitive_hash(
            "embedded interval backend hash",
            runtime.get("interval_backend_sha256"),
            PROFILE_INTERVAL_BACKEND,
        )
    elif dependency.name in {
        "ca_bridge_support_11_3299",
        "ca_direct_support_3299_56048351",
    }:
        configuration = report.get("configuration", {})
        coverage = report.get("profile_coverage", {})
        classification = report.get("classification", {})
        runtime = report.get("runtime_environment", {})
        expected_bounds = (
            (11, 3_299)
            if dependency.name == "ca_bridge_support_11_3299"
            else (3_299, 56_048_351)
        )
        record("certificate status", report.get("certificate_status") == "PASS", report.get("certificate_status"), "PASS")
        record("minimum support", configuration.get("min_support_prime") == expected_bounds[0], configuration.get("min_support_prime"), expected_bounds[0])
        record("maximum support", configuration.get("max_support_prime") == expected_bounds[1], configuration.get("max_support_prime"), expected_bounds[1])
        record("nonpositive profiles", coverage.get("nonpositive_or_unresolved_profile_count") == 0, coverage.get("nonpositive_or_unresolved_profile_count"), 0)
        record("classification status", classification.get("status") == "PASS", classification.get("status"), "PASS")
        record("all tie branches positive", classification.get("every_active_tie_has_two_positive_branches") is True, classification.get("every_active_tie_has_two_positive_branches"), True)
        embedded_hash = runtime.get("script_sha256")
    elif dependency.name == "ca_analytic_support_56048351_to_X":
        claim = report.get("claim", {})
        runtime = report.get("environment", {})
        named_checks = {
            row.get("name"): row for row in report.get("checks", [])
        }
        combined = named_checks.get(
            "minimum segmented combined clearance", {}
        )
        record("certificate status", report.get("status") == "PASS", report.get("status"), "PASS")
        expected_range = f"56048351 <= x <= {expected_x_max}"
        record("support range", claim.get("range") == expected_range, claim.get("range"), expected_range)
        record("profile scope", claim.get("profiles") == "every CA exponent profile with exact prime support x", claim.get("profiles"), "every CA exponent profile with exact prime support x")
        record("positive combined clearance", combined.get("status") == "PASS" and Decimal(str(combined.get("lower", "0"))) > Decimal("2e-7"), combined, "lower > 2e-7")
        embedded_hash = runtime.get("script_sha256")
        record_transitive_hash(
            "embedded residual interval backend hash",
            runtime.get("shared_interval_backend_sha256"),
            RESIDUAL_INTERVAL_BACKEND,
        )
        record_transitive_hash(
            "embedded previous upward certifier hash",
            runtime.get("previous_upward_certifier_sha256"),
            PREVIOUS_UPWARD_CERTIFIER,
        )
        dstar = report.get("independent_D_star_certificate", {})
        record_transitive_hash(
            "embedded D-star certifier hash",
            dstar.get("certifier_script_sha256"),
            DSTAR_CERTIFIER,
        )
        record_transitive_hash(
            "embedded D-star report hash",
            dstar.get("report_sha256"),
            DSTAR_REPORT,
        )
    else:
        raise ValueError(f"unknown dependency {dependency.name!r}")

    record("embedded source hash", embedded_hash == script_hash, embedded_hash, script_hash)
    summary = {
        "name": dependency.name,
        "status": "PASS" if all(row["status"] == "PASS" for row in checks) else "FAIL",
        "report_path": str(dependency.report.resolve()),
        "report_sha256": sha256(dependency.report),
        "script_path": str(dependency.script.resolve()),
        "script_sha256": script_hash,
        "transitive_files_checked": transitive_files,
        "checks": checks,
    }
    return report, summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate the certificate chain from the finite CA support theorem "
            "to an explicit all-integer Robin cutoff."
        )
    )
    parser.add_argument("--precision", type=int, default=60)
    parser.add_argument(
        "--x-max",
        type=int,
        default=None,
        help=(
            "Certified maximum CA support. By default it is parsed from the "
            "upward residual report and must match that report exactly."
        ),
    )
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser


def certified_x_max_from_upward_report() -> int:
    upward_path = DEPENDENCIES[-1].report
    report = json.loads(upward_path.read_bytes())
    range_text = str(report.get("claim", {}).get("range", ""))
    integers = [int(value) for value in re.findall(r"\d+", range_text)]
    if len(integers) != 2 or integers[0] != 56_048_351:
        raise ValueError(
            f"could not parse the certified support range from {range_text!r}"
        )
    return integers[1]


def clean_significant_digit_floor(value: int, significant_digits: int) -> int:
    """Round a positive integer down to a fixed number of significant digits."""

    if value <= 0:
        raise ValueError("the cutoff exponent must be positive")
    if significant_digits <= 0:
        raise ValueError("significant_digits must be positive")
    digits = len(str(value))
    if digits <= significant_digits:
        return value
    scale = 10 ** (digits - significant_digits)
    return value // scale * scale


def main() -> None:
    args = build_parser().parse_args()
    arithmetic = DecimalIntervals(args.precision)

    report_x_max = certified_x_max_from_upward_report()
    x_max = report_x_max if args.x_max is None else args.x_max
    if x_max != report_x_max:
        raise ValueError(
            f"x-max {x_max} does not match the upward report endpoint {report_x_max}"
        )
    y = x_max * SUPPORT_FRACTION_NUMERATOR // SUPPORT_FRACTION_DENOMINATOR
    if y < 3_275:
        raise ValueError("the 0.999X anchor lies below Dusart's x>=3275 domain")

    dependency_summaries: list[dict[str, Any]] = []
    for dependency in DEPENDENCIES:
        _, summary = load_dependency(dependency, x_max)
        dependency_summaries.append(summary)

    log_y = arithmetic.ln_integer(y)
    log_y_squared = arithmetic.mul(log_y, log_y)
    short_prime_relative_increment = arithmetic.div(
        arithmetic.point(1),
        arithmetic.mul(arithmetic.point(2), log_y_squared),
    )
    one_over_999 = arithmetic.rational(1, 999)
    short_increment_pass = short_prime_relative_increment[1] < one_over_999[0]
    short_prime_upper = arithmetic.mul(
        arithmetic.point(y),
        arithmetic.add(arithmetic.point(1), short_prime_relative_increment),
    )
    short_prime_below_x_pass = short_prime_upper[1] < Decimal(x_max)
    relative_theta_error = arithmetic.div(
        arithmetic.point(DUSART_THETA_CONSTANT), log_y_squared
    )
    absolute_theta_error = arithmetic.mul(
        arithmetic.point(y), relative_theta_error
    )
    theta_lower_expression = arithmetic.sub(
        arithmetic.point(y), absolute_theta_error
    )

    log_ten = arithmetic.ln_integer(10)
    directed_decimal_exponent = arithmetic.div(
        theta_lower_expression, log_ten
    )
    raw_floor_exponent = int(
        directed_decimal_exponent[0].to_integral_value(rounding=ROUND_FLOOR)
    )
    decimal_exponent = clean_significant_digit_floor(
        raw_floor_exponent, CLEAN_DECIMAL_EXPONENT_DIGITS
    )
    target_natural_log = arithmetic.mul(
        arithmetic.point(decimal_exponent), log_ten
    )
    cutoff_margin = arithmetic.floor.subtract(
        theta_lower_expression[0], target_natural_log[1]
    )
    numeric_cutoff_pass = (
        cutoff_margin > 0
        and short_increment_pass
        and short_prime_below_x_pass
    )
    dependencies_pass = all(
        row["status"] == "PASS" for row in dependency_summaries
    )
    overall_pass = dependencies_pass and numeric_cutoff_pass

    script_path = Path(__file__).resolve()
    report = {
        "certificate_status": "PASS" if overall_pass else "FAIL",
        "claim_certified_with_stated_external_theorems": (
            "Robin's inequality holds for every integer n with "
            f"5041 <= n <= 10^{decimal_exponent}."
        ),
        "ca_endpoint_definition": (
            "Dusart supplies a prime q in (floor(0.999X),X). Choose a "
            "generic CA parameter epsilon_0 strictly inside q's first-layer "
            "cell and outside the finite higher-layer threshold set. The "
            f"resulting unique CA profile C_q has exact support q<X={x_max}."
        ),
        "parametric_corollary": (
            "For any hash-validated certified support endpoint X in the same "
            "certificate chain, put y=floor(999X/1000). If y>=3275 and "
            "1/(2log(y)^2)<1/999, then there is a unique generic CA profile "
            "C_q of support q<X that exceeds "
            "exp(y*(1-3.965/log(y)^2)); a directed lower decimal exponent "
            "may be rounded down to state an ordinary-integer cutoff."
        ),
        "logical_chain": [
            "The directed base certificate proves Robin for every integer 5041<=n<=55440.",
            "The support-11 overlap certificate proves that the first exact-support-11 CA profile is 2^4*3^2*5*7*11=55440 (the old first-layer tie profile is 5040), so the base and CA chain meet at an identical endpoint.",
            f"The bridge, direct, and analytic CA certificates prove strict Robin for every CA profile with exact prime support 11<=P(C)<={x_max}, including all certified transition-tie branches.",
            "Put y=floor(0.999X). Dusart's short-prime interval gives a prime q with y<q<=y*(1+1/(2log(y)^2)). The directed check proves 1/(2log(y)^2)<1/999, so q<X.",
            "Choose epsilon_0 strictly inside q's first-layer cell and outside the finite set of higher-layer thresholds in that cell. The multiplicative CA rule then gives a unique profile C_q of exact support q.",
            "If epsilon_1>epsilon_2 and C_1,C_2 maximize sigma(n)/n^(1+epsilon_j), adding their two optimality inequalities gives (epsilon_1-epsilon_2)(log(C_2)-log(C_1))>=0. Hence CA size is nondecreasing as the parameter decreases. Since C_q is unique, every CA number below C_q belongs to a larger parameter and has support at most q<X.",
            "Robin 1984, Section 3, Proposition 1, gives f(n)<=max(f(C),f(C')) between consecutive CA numbers C<C', where f(n)=sigma(n)/(n log log n). Consequently the base and the complete preceding CA chain prove Robin through C_q.",
            "The profile C_q has full support, so log(C_q)>=theta(q)>theta(y).",
            f"Dusart 1999, Theorem 4, gives theta(y)>y*(1-3.965/log(y)^2) for y>=2. The directed numeric check divides this lower expression by log(10), rounds its integer lower bound down to {CLEAN_DECIMAL_EXPONENT_DIGITS} significant digits, and verifies a strict positive margin.",
            f"Therefore C_q>10^{decimal_exponent}, and the integer theorem follows.",
        ],
        "external_theorem_inputs_not_reproved_by_this_script": [
            {
                "theorem": "Robin consecutive-CA interpolation",
                "citation": "G. Robin, Grandes valeurs de la fonction somme des diviseurs et hypothese de Riemann, J. Math. Pures Appl. (9) 63 (1984), 187-213, Section 3, Proposition 1, p. 192, MR0774171.",
            },
            {
                "theorem": "Explicit theta estimate",
                "citation": "P. Dusart, Inegalites explicites pour psi(x), theta(x), pi(x) et les nombres premiers, C. R. Math. Rep. Acad. Sci. Canada 21 (1999), 53-59, Theorem 4.",
            },
            {
                "theorem": "Dusart short-prime interval",
                "citation": "P. Dusart, The kth prime is greater than k(ln k+ln ln k-1) for k>=2, Math. Comp. 68 (1999), 411-415, DOI 10.1090/S0025-5718-99-01037-6: for y>=3275 there is a prime y<p<=y(1+1/(2log(y)^2)).",
            },
            {
                "theorem": "CA threshold structure and full support",
                "citation": "L. Alaoglu and P. Erdos, On highly composite and similar numbers, Trans. Amer. Math. Soc. 56 (1944), 448-469; also recorded in Robin 1984, Section 2.",
            },
        ],
        "numeric_cutoff": {
            "X": str(x_max),
            "y_equals_floor_999X_over_1000": str(y),
            "support_fraction": "999/1000",
            "log_y_interval": interval_row(log_y),
            "short_prime_relative_increment_interval": interval_row(
                short_prime_relative_increment
            ),
            "required_short_prime_increment_upper_bound": "1/999",
            "short_prime_increment_check_status": (
                "PASS" if short_increment_pass else "FAIL"
            ),
            "short_prime_interval_upper_endpoint": interval_row(
                short_prime_upper
            ),
            "short_prime_upper_endpoint_is_below_X": (
                "PASS" if short_prime_below_x_pass else "FAIL"
            ),
            "dusart_relative_error_interval": interval_row(relative_theta_error),
            "theta_y_lower_expression_interval": interval_row(theta_lower_expression),
            "directed_decimal_exponent_interval": interval_row(
                directed_decimal_exponent
            ),
            "raw_floor_of_directed_lower_decimal_exponent": str(
                raw_floor_exponent
            ),
            "clean_significant_digits": CLEAN_DECIMAL_EXPONENT_DIGITS,
            "clean_decimal_exponent_K": str(
                decimal_exponent
            ),
            "K_log_10_interval": interval_row(target_natural_log),
            "strict_lower_margin_theta_expression_minus_K_log_10": str(cutoff_margin),
            "status": "PASS" if numeric_cutoff_pass else "FAIL",
        },
        "dependency_certificates": dependency_summaries,
        "rigor_boundary": [
            "This script validates the directed arithmetic and hashes the executable dependency chain; it does not reprove the cited external theorems.",
            "The upper endpoint is an ordinary-integer cutoff derived below the generic profile C_q. It is not the assertion that C_q itself equals a power of ten.",
            f"The conclusion is finite and does not extend the CA support theorem beyond X={x_max} or prove the Riemann hypothesis.",
        ],
        "runtime_environment": {
            "python_version": sys.version,
            "python_implementation": platform.python_implementation(),
            "decimal_module_version": decimal.__version__,
            "libmpdec_version": decimal.__libmpdec_version__,
            "precision_decimal_digits": args.precision,
            "script_path": str(script_path),
            "script_sha256": sha256(script_path),
            "argv": sys.argv,
        },
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(
        f"{report['certificate_status']}: dependencies="
        f"{'PASS' if dependencies_pass else 'FAIL'}, cutoff margin={cutoff_margin}; "
        f"report {args.report.resolve()}"
    )
    if not overall_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
