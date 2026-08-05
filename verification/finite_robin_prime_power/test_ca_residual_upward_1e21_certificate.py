from __future__ import annotations

import argparse
import decimal
import hashlib
import json
import platform
import sys
from decimal import Decimal
from pathlib import Path
from typing import Any

from test_ca_residual_finite_certificate import (
    DecimalIntervals,
    euler_gamma_interval,
    interval_row,
    machin_pi_interval,
    require_lower_above,
    require_upper_below,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_REPORT = ROOT / "ca_residual_upward_1e21_certificate_report.json"
BASE_CERTIFIER = ROOT / "test_ca_residual_finite_certificate.py"
D_STAR_REPORT = ROOT / "asymptotic_dstar_lower_bound_report.json"
D_STAR_CERTIFIER = ROOT / "test_asymptotic_dstar_lower_bound.py"

X_MIN = 56_048_351
X_SPLIT = 10**19
X_MAX = 10**21
ZERO_HEIGHT = 3_000_000_000_000
GLOBAL_H_COEFFICIENT = Decimal("2.035")
RESIDUAL_LOSS_TARGET = Decimal("0.457")
D_STAR_TARGET = Decimal("0.5")

Interval = tuple[Decimal, Decimal]


def load_asymptotic_dstar_certificate() -> dict[str, Any]:
    """Validate the independent analytic D_* tail certificate and source hash."""
    report_bytes = D_STAR_REPORT.read_bytes()
    report = json.loads(report_bytes)
    script_hash = hashlib.sha256(D_STAR_CERTIFIER.read_bytes()).hexdigest()
    embedded_hash = report["runtime_environment"]["certifier_script_sha256"]
    if report["status"] != "PASS":
        raise AssertionError("the asymptotic D_* certificate is not PASS")
    if embedded_hash != script_hash:
        raise AssertionError(
            "the asymptotic D_* report does not match its current certifier script"
        )
    if report["endpoint"]["x0"] != X_MIN:
        raise AssertionError("the asymptotic D_* certificate starts at the wrong endpoint")
    if Decimal(report["parameters"]["target"]) != D_STAR_TARGET:
        raise AssertionError("the asymptotic D_* certificate used a different target")
    if not all(report["monotonicity_checks"].values()):
        raise AssertionError("an asymptotic D_* monotonicity check failed")
    minimum_interval = (
        Decimal(report["endpoint"]["scaled_dstar_lower"]["lower"]),
        Decimal(report["endpoint"]["scaled_dstar_lower"]["upper"]),
    )
    if minimum_interval[0] <= D_STAR_TARGET:
        raise AssertionError("the asymptotic D_* lower bound is not above 0.5")
    return {
        "report": report,
        "minimum_interval": minimum_interval,
        "script_sha256": script_hash,
        "report_sha256": hashlib.sha256(report_bytes).hexdigest(),
    }


def q_interval(arithmetic: DecimalIntervals, x: int) -> Interval:
    """Return Q(x)=x*(1+1/(2*log(x)^2))+1 with outward rounding."""
    one = arithmetic.point(1)
    two = arithmetic.point(2)
    x_interval = arithmetic.point(x)
    log_x = arithmetic.ln(x_interval)
    correction = arithmetic.div(one, arithmetic.mul(two, arithmetic.mul(log_x, log_x)))
    return arithmetic.add(arithmetic.mul(x_interval, arithmetic.add(one, correction)), one)


def layer_majorant(
    arithmetic: DecimalIntervals, q_value: Interval
) -> tuple[Interval, Interval]:
    """Return the analytic majorant for F(Q)/sqrt(Q) and its tail derivative margin."""
    one = arithmetic.point(1)
    two = arithmetic.point(2)
    three = arithmetic.point(3)
    four = arithmetic.point(4)
    five = arithmetic.point(5)
    log_two = arithmetic.ln(two)
    log_q = arithmetic.ln(q_value)
    log2_q = arithmetic.div(log_q, log_two)
    b_q = arithmetic.add(
        log2_q,
        arithmetic.div(arithmetic.ln(arithmetic.mul(two, log2_q)), log_two),
    )
    b_minus_four = arithmetic.sub(b_q, four)
    majorant = arithmetic.add(
        arithmetic.add(
            arithmetic.sqrt(two),
            arithmetic.mul(
                arithmetic.rational_power(three, 1, 3),
                arithmetic.rational_power(q_value, -1, 6),
            ),
        ),
        arithmetic.add(
            arithmetic.mul(
                arithmetic.rational_power(four, 1, 4),
                arithmetic.rational_power(q_value, -1, 4),
            ),
            arithmetic.mul(
                arithmetic.mul(arithmetic.rational_power(five, 1, 5), b_minus_four),
                arithmetic.rational_power(q_value, -3, 10),
            ),
        ),
    )
    b_derivative = arithmetic.div(
        arithmetic.add(one, arithmetic.div(one, log_q)), log_two
    )
    derivative_margin = arithmetic.sub(
        arithmetic.mul(arithmetic.rational(3, 10), b_minus_four), b_derivative
    )
    return majorant, derivative_margin


def partial_rh_z_data(
    arithmetic: DecimalIntervals, x: int, pi: Interval
) -> dict[str, Interval]:
    """Compute the Buethe-2016 theta coefficient and the resulting z enclosure."""
    one = arithmetic.point(1)
    two = arithmetic.point(2)
    four = arithmetic.point(4)
    eight = arithmetic.point(8)
    x_value = arithmetic.point(x)
    log_x = arithmetic.ln(x_value)
    sqrt_x = arithmetic.sqrt(x_value)
    theta_coefficient = arithmetic.div(
        arithmetic.mul(log_x, arithmetic.sub(log_x, two)),
        arithmetic.mul(eight, pi),
    )
    absolute_z_coefficient = arithmetic.add(
        theta_coefficient, arithmetic.point(GLOBAL_H_COEFFICIENT)
    )
    delta = arithmetic.div(absolute_z_coefficient, sqrt_x)
    one_minus_delta = arithmetic.sub(one, delta)
    a_minus = arithmetic.add(log_x, arithmetic.ln(one_minus_delta))
    curvature = arithmetic.div(
        arithmetic.add(a_minus, one),
        arithmetic.mul(
            arithmetic.mul(one_minus_delta, one_minus_delta),
            arithmetic.mul(a_minus, a_minus),
        ),
    )
    b2_scaled_loss = arithmetic.mul(
        arithmetic.mul(sqrt_x, log_x),
        arithmetic.mul(
            arithmetic.div(arithmetic.mul(delta, delta), two), curvature
        ),
    )
    height_requirement = arithmetic.mul(
        arithmetic.point("4.92"), arithmetic.sqrt(arithmetic.div(x_value, log_x))
    )
    return {
        "log_x": log_x,
        "sqrt_x": sqrt_x,
        "theta_coefficient": theta_coefficient,
        "absolute_z_coefficient": absolute_z_coefficient,
        "delta": delta,
        "a_minus": a_minus,
        "b2_scaled_loss": b2_scaled_loss,
        "height_requirement": height_requirement,
    }


def low_zero_loss(
    arithmetic: DecimalIntervals,
    log_x: Interval,
    sqrt_x: Interval,
    zero_sum: Interval,
    log_two_pi: Interval,
) -> Interval:
    one = arithmetic.point(1)
    three = arithmetic.point(3)
    four = arithmetic.point(4)
    kernel_factor = arithmetic.add(
        one,
        arithmetic.add(
            arithmetic.div(three, log_x),
            arithmetic.div(four, arithmetic.mul(log_x, log_x)),
        ),
    )
    return arithmetic.add(
        arithmetic.mul(zero_sum, kernel_factor),
        arithmetic.div(log_two_pi, sqrt_x),
    )


def high_zero_loss(
    arithmetic: DecimalIntervals,
    log_x: Interval,
    sqrt_x: Interval,
    log_zero_height: Interval,
    pi: Interval,
) -> Interval:
    one = arithmetic.point(1)
    four = arithmetic.point(4)
    zero_height = arithmetic.point(ZERO_HEIGHT)
    return arithmetic.div(
        arithmetic.mul(
            arithmetic.mul(
                arithmetic.mul(four, sqrt_x),
                arithmetic.div(arithmetic.add(log_x, one), log_x),
            ),
            arithmetic.add(log_zero_height, one),
        ),
        arithmetic.mul(pi, zero_height),
    )


def build_certificate(precision: int, gamma_terms: int) -> dict[str, Any]:
    arithmetic = DecimalIntervals(precision)
    dstar_certificate = load_asymptotic_dstar_certificate()
    dstar_report = dstar_certificate["report"]
    dstar_minimum = dstar_certificate["minimum_interval"]
    one = arithmetic.point(1)
    two = arithmetic.point(2)
    four = arithmetic.point(4)

    pi = machin_pi_interval(arithmetic)
    gamma = euler_gamma_interval(arithmetic, gamma_terms)
    zero_sum = arithmetic.sub(
        arithmetic.add(two, gamma), arithmetic.ln(arithmetic.mul(four, pi))
    )
    log_two_pi = arithmetic.ln(arithmetic.mul(two, pi))

    q_min = q_interval(arithmetic, X_MIN)
    q_max = q_interval(arithmetic, X_MAX)
    q_over_x_min = arithmetic.div(q_min, arithmetic.point(X_MIN))
    normalized_layer_majorant, layer_tail_derivative_margin = layer_majorant(
        arithmetic, q_min
    )
    twice_sqrt_q_over_x = arithmetic.mul(two, arithmetic.sqrt(q_over_x_min))
    rosser_h_coefficient = arithmetic.mul(
        arithmetic.point("1.01624"), twice_sqrt_q_over_x
    )
    finite_theta_h_coefficient = twice_sqrt_q_over_x
    maximum_layer_cutoff = arithmetic.sqrt(arithmetic.mul(two, q_max))

    left = partial_rh_z_data(arithmetic, X_MIN, pi)
    split = partial_rh_z_data(arithmetic, X_SPLIT, pi)
    right = partial_rh_z_data(arithmetic, X_MAX, pi)

    low_left = low_zero_loss(
        arithmetic, left["log_x"], left["sqrt_x"], zero_sum, log_two_pi
    )
    low_split = low_zero_loss(
        arithmetic, split["log_x"], split["sqrt_x"], zero_sum, log_two_pi
    )

    zero_height = arithmetic.point(ZERO_HEIGHT)
    log_zero_height = arithmetic.ln(zero_height)
    e_interval = arithmetic.exp(one)
    log_two_pi_e = arithmetic.ln(arithmetic.mul(arithmetic.mul(two, pi), e_interval))
    zero_count_error = arithmetic.add(
        arithmetic.add(
            arithmetic.mul(arithmetic.point("0.1038"), log_zero_height),
            arithmetic.mul(arithmetic.point("0.2573"), arithmetic.ln(log_zero_height)),
        ),
        arithmetic.point("9.3675"),
    )
    zero_count_margin_at_t = arithmetic.sub(
        arithmetic.div(arithmetic.mul(zero_height, log_two_pi_e), arithmetic.mul(two, pi)),
        zero_count_error,
    )
    zero_count_derivative_margin = arithmetic.sub(
        arithmetic.div(log_two_pi_e, arithmetic.mul(two, pi)),
        arithmetic.add(
            arithmetic.div(arithmetic.point("0.1038"), zero_height),
            arithmetic.div(
                arithmetic.point("0.2573"), arithmetic.mul(zero_height, log_zero_height)
            ),
        ),
    )
    high_split = high_zero_loss(
        arithmetic, split["log_x"], split["sqrt_x"], log_zero_height, pi
    )
    high_right = high_zero_loss(
        arithmetic, right["log_x"], right["sqrt_x"], log_zero_height, pi
    )

    # On each segment, combine the decreasing low/B2 losses at the left endpoint
    # with the increasing high-zero loss at the right endpoint.
    lower_segment_loss = arithmetic.add(
        arithmetic.add(low_left, left["b2_scaled_loss"]), high_split
    )
    upper_segment_loss = arithmetic.add(
        arithmetic.add(low_split, split["b2_scaled_loss"]), high_right
    )
    uniform_residual_loss = (
        Decimal(0), max(lower_segment_loss[1], upper_segment_loss[1])
    )
    conservative_clearance = arithmetic.sub(
        arithmetic.point(D_STAR_TARGET), uniform_residual_loss
    )
    certified_clearance = arithmetic.sub(dstar_minimum, uniform_residual_loss)

    # The derivative proof for the B2 endpoint reduction uses
    # 1/L+2c'/c <= (5L-6)/(L(L-2)) < 1/2 at the left endpoint.
    b2_log_derivative_upper = arithmetic.div(
        arithmetic.sub(arithmetic.mul(arithmetic.point(5), left["log_x"]), arithmetic.point(6)),
        arithmetic.mul(left["log_x"], arithmetic.sub(left["log_x"], two)),
    )
    b2_monotonicity_margin = arithmetic.sub(
        arithmetic.rational(1, 2), b2_log_derivative_upper
    )
    delta_log_derivative_upper = arithmetic.div(
        arithmetic.mul(two, arithmetic.sub(left["log_x"], one)),
        arithmetic.mul(left["log_x"], arithmetic.sub(left["log_x"], two)),
    )
    delta_monotonicity_margin = arithmetic.sub(
        arithmetic.rational(1, 2), delta_log_derivative_upper
    )

    domain_checks = {
        "Dusart_x_min_at_least_3275": X_MIN >= 3_275,
        "Buethe_2016_equation_7_4_x_min_at_least_5000": X_MIN >= 5_000,
        "Buethe_2018_layer_cutoff_domain_x_at_most_1e19": maximum_layer_cutoff[1]
        < Decimal("1e19"),
        "Platt_Trudgian_height_covers_Buethe_requirement": right["height_requirement"][1]
        < Decimal(ZERO_HEIGHT),
        "HSW_zero_count_height_at_least_e": Decimal(ZERO_HEIGHT) > e_interval[1],
    }
    if not all(domain_checks.values()):
        raise AssertionError(f"domain checks failed: {domain_checks}")

    checks = [
        require_upper_below(
            "normalized analytic layer majorant at X_MIN",
            normalized_layer_majorant,
            Decimal(2),
        ),
        require_lower_above(
            "layer-tail monotonicity margin",
            layer_tail_derivative_margin,
            Decimal(0),
        ),
        require_upper_below(
            "global Rosser-Schoenfeld H coefficient",
            rosser_h_coefficient,
            GLOBAL_H_COEFFICIENT,
        ),
        require_upper_below(
            "stronger finite-theta H coefficient",
            finite_theta_h_coefficient,
            Decimal("2.002"),
        ),
        require_upper_below(
            "largest layer cutoff", maximum_layer_cutoff, Decimal("1e19")
        ),
        require_upper_below(
            "partial-RH height requirement at X_MAX",
            right["height_requirement"],
            Decimal(ZERO_HEIGHT),
        ),
        require_lower_above(
            "delta monotonicity margin", delta_monotonicity_margin, Decimal(0)
        ),
        require_lower_above(
            "B2 monotonicity margin", b2_monotonicity_margin, Decimal(0)
        ),
        require_upper_below(
            "delta at X_MIN is below one", left["delta"], Decimal(1)
        ),
        require_lower_above(
            "A_minus at X_MIN is positive", left["a_minus"], Decimal(0)
        ),
        require_lower_above(
            "zero-count upper-bound margin at T", zero_count_margin_at_t, Decimal(0)
        ),
        require_lower_above(
            "zero-count derivative margin", zero_count_derivative_margin, Decimal(0)
        ),
        require_upper_below(
            "lower-segment residual loss", lower_segment_loss, Decimal("0.109")
        ),
        require_upper_below(
            "upper-segment residual loss", upper_segment_loss, RESIDUAL_LOSS_TARGET
        ),
        require_upper_below(
            "uniform upward residual loss", uniform_residual_loss, RESIDUAL_LOSS_TARGET
        ),
        require_lower_above(
            "conservative combined clearance from D_* target 0.5",
            conservative_clearance,
            Decimal("0.043"),
        ),
        require_lower_above(
            "combined clearance from certified D_* lower bound",
            certified_clearance,
            Decimal("0.059"),
        ),
    ]

    script_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    base_hash = hashlib.sha256(BASE_CERTIFIER.read_bytes()).hexdigest()
    return {
        "status": "PASS",
        "claim": {
            "range": f"{X_MIN} <= x <= {X_MAX}",
            "profiles": "every CA exponent profile with exact prime support x",
            "residual_bound": (
                "sqrt(x)log(x)*(R_infinity(x)+B2(n,x))>-0.457"
            ),
            "combined_finite_CA_theorem": (
                "the hash-validated asymptotic D_* certificate proves "
                "sqrt(x)log(x)D_*(x)>0.5, hence sqrt(x)log(x)G(n,x)>0.043"
            ),
        },
        "independent_D_star_certificate": {
            "status": dstar_report["status"],
            "claim": dstar_report["claim"],
            "report_path": str(D_STAR_REPORT),
            "report_sha256": dstar_certificate["report_sha256"],
            "certifier_path": str(D_STAR_CERTIFIER),
            "certifier_script_sha256": dstar_certificate["script_sha256"],
            "embedded_certifier_script_sha256": dstar_report["runtime_environment"][
                "certifier_script_sha256"
            ],
            "endpoint": dstar_report["endpoint"]["x0"],
            "uniform_scaled_D_star_lower": interval_row(dstar_minimum),
            "strict_target": dstar_report["parameters"]["target"],
            "method": dstar_report["method"],
        },
        "external_theorem_inputs": [
            {
                "source": "Dusart (1999)",
                "input": "q<=x(1+1/(2log(x)^2)) for x>=3275",
            },
            {
                "source": "Rosser-Schoenfeld (1962), equation (3.24)",
                "input": "theta(t)<1.01624t for every t>0",
            },
            {
                "source": "Buethe, Math. Comp. 85 (2016), Theorem 2 proof equation (7.4)",
                "input": (
                    "under partial RH through T and 4.92sqrt(x/log x)<=T, "
                    "|theta(x)-x|<=sqrt(x)log(x)(log(x)-2)/(8pi) for x>=5000"
                ),
            },
            {
                "source": "Buethe, Math. Comp. 87 (2018), Theorem 2",
                "input": "theta(u)<u for 1<=u<=1e19; used only for the stronger layer audit",
            },
            {
                "source": "Platt-Trudgian (2021)",
                "input": "all zeta zeros through height 3e12 satisfy RH",
            },
            {
                "source": "Hasanalizade-Shen-Wong (2022), Corollary 1.2",
                "input": "explicit zero-count estimate for t>=e",
            },
        ],
        "analytic_lemmas": {
            "layer_envelope": (
                "The analytic F(Q)/sqrt(Q) majorant decreases from X_MIN and is <2. "
                "Also Q(x)/x decreases, giving H<2.035sqrt(x) globally."
            ),
            "maximum_layer_cutoff": (
                "For Q>e/2, (sQ)^(1/s) decreases for s>=2; the largest cutoff is sqrt(2Q)."
            ),
            "partial_RH_domain": (
                "sqrt(x/log x) increases for log x>1, so the height requirement is maximal at X_MAX."
            ),
            "b2_monotonicity": (
                "Let c(L)=L(L-2)/(8pi)+2.035 and delta=c(L)e^(-L/2). "
                "Both delta and the curvature factor decrease. Moreover "
                "1/L+2c'/c<1/2, so e^(-L/2)Lc(L)^2 decreases."
            ),
            "low_high_monotonicity": (
                "The low-zero loss decreases; sqrt(x)(L+1)/L and hence the high-zero loss increase."
            ),
            "termwise_zero_integration": (
                "Symmetric finite zero truncations are integrated first. The integrated "
                "tail is absolutely convergent because |J_rho/rho|<=2h(L)/gamma^2 "
                "and N(t)=O(tlog t)."
            ),
            "two_segment_cover": (
                f"On [{X_MIN},{X_SPLIT}] use low/B2 at the left and high at the right; "
                f"on [{X_SPLIT},{X_MAX}] do the same."
            ),
        },
        "constants": {
            "pi": interval_row(pi),
            "euler_gamma": interval_row(gamma),
            "exact_zero_sum": interval_row(zero_sum),
            "Q_at_X_MIN": interval_row(q_min),
            "Q_at_X_MAX": interval_row(q_max),
            "normalized_layer_majorant": interval_row(normalized_layer_majorant),
            "layer_tail_derivative_margin": interval_row(layer_tail_derivative_margin),
            "global_H_coefficient": interval_row(rosser_h_coefficient),
            "finite_theta_H_coefficient": interval_row(finite_theta_h_coefficient),
            "maximum_layer_cutoff": interval_row(maximum_layer_cutoff),
            "partial_RH_height_requirement_at_X_MAX": interval_row(
                right["height_requirement"]
            ),
            "B2_monotonicity_margin": interval_row(b2_monotonicity_margin),
            "delta_monotonicity_margin": interval_row(delta_monotonicity_margin),
            "B2_loss_at_X_MIN": interval_row(left["b2_scaled_loss"]),
            "B2_loss_at_split": interval_row(split["b2_scaled_loss"]),
            "low_zero_loss_at_X_MIN": interval_row(low_left),
            "low_zero_loss_at_split": interval_row(low_split),
            "high_zero_loss_at_split": interval_row(high_split),
            "high_zero_loss_at_X_MAX": interval_row(high_right),
            "lower_segment_loss": interval_row(lower_segment_loss),
            "upper_segment_loss": interval_row(upper_segment_loss),
            "uniform_residual_loss": interval_row(uniform_residual_loss),
            "conservative_clearance_from_D_star_target_0_5": interval_row(
                conservative_clearance
            ),
            "clearance_from_certified_D_star_lower": interval_row(certified_clearance),
            "zero_count_margin_at_T": interval_row(zero_count_margin_at_t),
            "zero_count_derivative_margin": interval_row(zero_count_derivative_margin),
        },
        "checks": checks,
        "exact_domain_checks": domain_checks,
        "rigor_scope": {
            "proved_here": (
                "All layer, domain, monotonic endpoint, two-segment, and scalar "
                "comparisons recorded in this report."
            ),
            "external": (
                "The six published theorem inputs are accepted and not recomputed."
            ),
            "combined_status": (
                "The independent asymptotic D_* certificate passed and its embedded "
                "source hash matches, so the combined CA theorem is complete on the range."
            ),
            "remaining_gap": (
                "This is not a proof-assistant formalization, and it does not perform "
                "the reduction from CA supports to every integer. The published inputs "
                "and both executable certificate runtimes remain trusted."
            ),
            "decimal_assumption": (
                "CPython Decimal transcendental operations have documented "
                "correct-rounding semantics; the shared backend pads every such endpoint."
            ),
        },
        "environment": {
            "python": sys.version,
            "implementation": platform.python_implementation(),
            "decimal_module": decimal.__file__,
            "libmpdec": getattr(decimal, "__libmpdec_version__", None),
            "precision": precision,
            "gamma_harmonic_terms": gamma_terms,
            "script_sha256": script_hash,
            "shared_interval_backend_path": str(BASE_CERTIFIER),
            "shared_interval_backend_sha256": base_hash,
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Certify the CA residual bound from 56,048,351 through 1e21."
    )
    parser.add_argument("--precision", type=int, default=60)
    parser.add_argument("--gamma-terms", type=int, default=100_000)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    certificate = build_certificate(args.precision, args.gamma_terms)
    args.report.write_text(json.dumps(certificate, indent=2), encoding="utf-8")
    constants = certificate["constants"]
    print(f"status: {certificate['status']}")
    print(f"global H coefficient upper: {constants['global_H_coefficient']['upper']}")
    print(
        "partial-RH requirement upper: "
        f"{constants['partial_RH_height_requirement_at_X_MAX']['upper']}"
    )
    print(f"B2 loss at X_MIN upper: {constants['B2_loss_at_X_MIN']['upper']}")
    print(f"high-zero loss at X_MAX upper: {constants['high_zero_loss_at_X_MAX']['upper']}")
    print(f"lower-segment loss upper: {constants['lower_segment_loss']['upper']}")
    print(f"upper-segment loss upper: {constants['upper_segment_loss']['upper']}")
    print(f"uniform residual loss upper: {constants['uniform_residual_loss']['upper']}")
    print(
        "conservative combined clearance lower: "
        f"{constants['conservative_clearance_from_D_star_target_0_5']['lower']}"
    )
    print(
        "certified-bound combined clearance lower: "
        f"{constants['clearance_from_certified_D_star_lower']['lower']}"
    )
    print(f"wrote {args.report}")


if __name__ == "__main__":
    main()
