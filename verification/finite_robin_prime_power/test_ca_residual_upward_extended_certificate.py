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
from test_ca_residual_upward_1e21_certificate import (
    D_STAR_CERTIFIER,
    D_STAR_REPORT,
    GLOBAL_H_COEFFICIENT,
    layer_majorant,
    load_asymptotic_dstar_certificate,
    low_zero_loss,
    partial_rh_z_data,
    q_interval,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_REPORT = ROOT / "ca_residual_upward_extended_certificate_report.json"
BASE_CERTIFIER = ROOT / "test_ca_residual_finite_certificate.py"
PREVIOUS_UPWARD_CERTIFIER = ROOT / "test_ca_residual_upward_1e21_certificate.py"

X_MIN = 56_048_351
X_MAX = 164_967_000_000_000_000_000_000
X_FINE_START = 100_000_000_000_000_000_000_000
X_FINE_STEP = 100_000_000_000_000_000_000
ZERO_HEIGHT = 3_000_000_000_000
DUSART_ETA = Decimal("3.965")
SEGMENT_CLEARANCE_TARGET = Decimal("0.0000002")

# These two points rigorously bracket the sign change of the scalar pointwise
# envelope used by this certificate.  A negative value at the upper point is
# not a counterexample; it only records the limit of this particular envelope.
BARRIER_LOW = 164_967_000_000_000_000_000_000
BARRIER_HIGH = 164_968_000_000_000_000_000_000

HSW_A = Decimal("0.1038")
HSW_B = Decimal("0.2573")
HSW_C = Decimal("9.3675")

Interval = tuple[Decimal, Decimal]


def dstar_closed_lower(arithmetic: DecimalIntervals, x: int) -> dict[str, Interval]:
    """Evaluate the increasing closed lower bound for sqrt(x) log(x) D_*(x)."""
    one = arithmetic.point(1)
    two = arithmetic.point(2)
    three = arithmetic.point(3)
    x_value = arithmetic.point(x)
    log_x = arithmetic.ln(x_value)
    log_two = arithmetic.ln(two)
    ell = arithmetic.div(arithmetic.add(log_x, log_two), two)
    sqrt_two = arithmetic.sqrt(two)
    sqrt_x = arithmetic.sqrt(x_value)
    y = arithmetic.sqrt(arithmetic.mul(two, x_value))
    r = arithmetic.sqrt(arithmetic.div(x_value, two))

    epsilon = arithmetic.div(
        arithmetic.point(DUSART_ETA), arithmetic.mul(ell, ell)
    )
    a_term = arithmetic.div(log_x, arithmetic.mul(sqrt_two, ell))
    q = arithmetic.sub(one, arithmetic.div(one, r))
    j_lower = arithmetic.mul(
        arithmetic.div(log_x, sqrt_two),
        arithmetic.sub(
            arithmetic.div(q, ell),
            arithmetic.div(one, arithmetic.mul(ell, ell)),
        ),
    )
    scaled_s2_lower = arithmetic.sub(
        arithmetic.sub(
            arithmetic.mul(arithmetic.sub(one, epsilon), j_lower),
            arithmetic.mul(two, arithmetic.mul(epsilon, a_term)),
        ),
        arithmetic.div(one, sqrt_x),
    )
    alpha = arithmetic.sub(
        arithmetic.rational(1, 2), arithmetic.div(one, arithmetic.mul(three, y))
    )
    scaled_dstar_lower = arithmetic.mul(alpha, scaled_s2_lower)
    return {
        "log_x": log_x,
        "ell": ell,
        "sqrt_x": sqrt_x,
        "y": y,
        "r": r,
        "epsilon": epsilon,
        "A": a_term,
        "q": q,
        "J_lower": j_lower,
        "scaled_prime_square_sum_lower": scaled_s2_lower,
        "alpha": alpha,
        "scaled_dstar_lower": scaled_dstar_lower,
    }


def hsw_tail_sums(
    arithmetic: DecimalIntervals, pi: Interval
) -> dict[str, Interval]:
    """Bound the positive-ordinate sums S_2(T) and S_3(T) from HSW."""
    one = arithmetic.point(1)
    two = arithmetic.point(2)
    zero_height = arithmetic.point(ZERO_HEIGHT)
    log_t = arithmetic.ln(zero_height)
    log_log_t = arithmetic.ln(log_t)

    # For M(t)=t/(2*pi) log(t/(2*pi*e)) and |N(t)-M(t)|<=E(t),
    # Stieltjes integration gives
    # S_2=-N(T)/T^2+2 integral_T^infinity N(t)/t^3 dt.
    # The main term is exact.  The log-log error integral uses
    # 1/log(t)<=1/log(T), which is valid throughout the tail.
    main = arithmetic.div(
        arithmetic.add(
            arithmetic.ln(
                arithmetic.div(zero_height, arithmetic.mul(two, pi))
            ),
            one,
        ),
        arithmetic.mul(arithmetic.mul(two, pi), zero_height),
    )
    error_numerator = arithmetic.add(
        arithmetic.add(
            arithmetic.mul(
                arithmetic.mul(two, arithmetic.point(HSW_A)), log_t
            ),
            arithmetic.div(arithmetic.point(HSW_A), two),
        ),
        arithmetic.add(
            arithmetic.add(
                arithmetic.mul(
                    arithmetic.mul(two, arithmetic.point(HSW_B)), log_log_t
                ),
                arithmetic.div(
                    arithmetic.point(HSW_B), arithmetic.mul(two, log_t)
                ),
            ),
            arithmetic.mul(two, arithmetic.point(HSW_C)),
        ),
    )
    error = arithmetic.div(
        error_numerator, arithmetic.mul(zero_height, zero_height)
    )
    s2_upper = arithmetic.add(main, error)
    s3_upper = arithmetic.div(s2_upper, zero_height)
    return {
        "log_T": log_t,
        "main": main,
        "error": error,
        "S2_upper": s2_upper,
        "S3_upper": s3_upper,
    }


def sharpened_high_zero_loss(
    arithmetic: DecimalIntervals,
    x: int,
    tail_sums: dict[str, Interval],
) -> Interval:
    """Return the scaled high-zero loss after two integrations by parts."""
    one = arithmetic.point(1)
    two = arithmetic.point(2)
    x_value = arithmetic.point(x)
    log_x = arithmetic.ln(x_value)
    sqrt_x = arithmetic.sqrt(x_value)
    log_x_squared = arithmetic.mul(log_x, log_x)
    h = arithmetic.div(arithmetic.add(log_x, one), log_x_squared)
    minus_h_prime = arithmetic.add(
        arithmetic.div(one, log_x_squared),
        arithmetic.div(two, arithmetic.mul(log_x_squared, log_x)),
    )
    kernel_tail = arithmetic.add(
        arithmetic.mul(h, tail_sums["S2_upper"]),
        arithmetic.mul(
            arithmetic.mul(two, minus_h_prime), tail_sums["S3_upper"]
        ),
    )
    reflection_factor = arithmetic.add(one, arithmetic.div(one, x_value))
    return arithmetic.mul(
        arithmetic.mul(arithmetic.mul(sqrt_x, log_x), reflection_factor),
        kernel_tail,
    )


def pointwise_components(
    arithmetic: DecimalIntervals,
    x: int,
    pi: Interval,
    zero_sum: Interval,
    log_two_pi: Interval,
    tail_sums: dict[str, Interval],
) -> dict[str, Interval]:
    z_data = partial_rh_z_data(arithmetic, x, pi)
    dstar_data = dstar_closed_lower(arithmetic, x)
    low_loss = low_zero_loss(
        arithmetic,
        z_data["log_x"],
        z_data["sqrt_x"],
        zero_sum,
        log_two_pi,
    )
    high_loss = sharpened_high_zero_loss(arithmetic, x, tail_sums)
    total_loss = arithmetic.add(
        arithmetic.add(low_loss, z_data["b2_scaled_loss"]), high_loss
    )
    clearance = arithmetic.sub(dstar_data["scaled_dstar_lower"], total_loss)
    return {
        "dstar": dstar_data["scaled_dstar_lower"],
        "low_zero_loss": low_loss,
        "b2_loss": z_data["b2_scaled_loss"],
        "high_zero_loss": high_loss,
        "total_loss": total_loss,
        "clearance": clearance,
        "height_requirement": z_data["height_requirement"],
    }


def segment_boundaries() -> list[int]:
    coarse = [
        X_MIN,
        10**9,
        10**11,
        10**13,
        10**15,
        10**17,
        10**19,
        10**20,
        10**21,
        10**22,
        5 * 10**22,
        X_FINE_START,
    ]
    fine = list(range(X_FINE_START + X_FINE_STEP, X_MAX, X_FINE_STEP))
    # The tiny final segment preserves a directed clearance near the scalar
    # envelope barrier without making the complete partition unnecessarily fine.
    boundaries = coarse + fine + [
        164_950_000_000_000_000_000_000,
        X_MAX,
    ]
    if boundaries[-1] != X_MAX or any(
        left >= right for left, right in zip(boundaries, boundaries[1:])
    ):
        raise AssertionError("invalid segment boundary construction")
    return boundaries


def build_certificate(precision: int, gamma_terms: int) -> dict[str, Any]:
    arithmetic = DecimalIntervals(precision)
    one = arithmetic.point(1)
    two = arithmetic.point(2)
    four = arithmetic.point(4)
    pi = machin_pi_interval(arithmetic)
    gamma = euler_gamma_interval(arithmetic, gamma_terms)
    zero_sum = arithmetic.sub(
        arithmetic.add(two, gamma), arithmetic.ln(arithmetic.mul(four, pi))
    )
    log_two_pi = arithmetic.ln(arithmetic.mul(two, pi))
    dstar_certificate = load_asymptotic_dstar_certificate()
    dstar_report = dstar_certificate["report"]
    tail_sums = hsw_tail_sums(arithmetic, pi)

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
    maximum_layer_cutoff = arithmetic.sqrt(arithmetic.mul(two, q_max))

    # High-zero monotonicity: after scaling, the bound is
    # (e^(L/2)+e^(-L/2))*A(L), where the factor L is already absorbed:
    # A(L)=S2*(1+1/L)+S3*(2/L+4/L^2).  Thus -A'/A (not -h'/h) is bounded by
    # 1/L^2+T^-1*(2/L^2+8/L^3), since A>=S2 and S3<=S2/T.  The base term
    # grows and this rational upper bound decreases.
    dstar_at_min = dstar_closed_lower(arithmetic, X_MIN)
    log_x_min = dstar_at_min["log_x"]
    left_z = partial_rh_z_data(arithmetic, X_MIN, pi)
    inv_x_min = arithmetic.div(one, arithmetic.point(X_MIN))
    high_base_derivative = arithmetic.mul(
        arithmetic.rational(1, 2),
        arithmetic.div(
            arithmetic.sub(one, inv_x_min), arithmetic.add(one, inv_x_min)
        ),
    )
    log_x_min_squared = arithmetic.mul(log_x_min, log_x_min)
    high_amplitude_derivative = arithmetic.add(
        arithmetic.div(one, log_x_min_squared),
        arithmetic.div(
            arithmetic.add(
                arithmetic.div(two, log_x_min_squared),
                arithmetic.div(
                    arithmetic.point(8),
                    arithmetic.mul(log_x_min_squared, log_x_min),
                ),
            ),
            arithmetic.point(ZERO_HEIGHT),
        ),
    )
    high_monotonicity_margin = arithmetic.sub(
        high_base_derivative, high_amplitude_derivative
    )

    # Recheck the two endpoint reductions used for B2 instead of relying only
    # on the earlier report.  For c(L)=L(L-2)/(8*pi)+2.035, these rational
    # upper bounds control respectively delta'/delta and the remaining
    # normalized B2 factor.
    delta_log_derivative_upper = arithmetic.div(
        arithmetic.mul(two, arithmetic.sub(log_x_min, one)),
        arithmetic.mul(log_x_min, arithmetic.sub(log_x_min, two)),
    )
    delta_monotonicity_margin = arithmetic.sub(
        arithmetic.rational(1, 2), delta_log_derivative_upper
    )
    b2_log_derivative_upper = arithmetic.div(
        arithmetic.sub(
            arithmetic.mul(arithmetic.point(5), log_x_min), arithmetic.point(6)
        ),
        arithmetic.mul(log_x_min, arithmetic.sub(log_x_min, two)),
    )
    b2_monotonicity_margin = arithmetic.sub(
        arithmetic.rational(1, 2), b2_log_derivative_upper
    )

    # The stronger variable-q D_* expression agrees with the independent
    # endpoint report at X_MIN and is increasing by the same elementary signs.
    reported_dstar_min = (
        Decimal(dstar_report["endpoint"]["scaled_dstar_lower"]["lower"]),
        Decimal(dstar_report["endpoint"]["scaled_dstar_lower"]["upper"]),
    )
    dstar_endpoint_overlap = not (
        dstar_at_min["scaled_dstar_lower"][1] < reported_dstar_min[0]
        or reported_dstar_min[1] < dstar_at_min["scaled_dstar_lower"][0]
    )
    log_two = arithmetic.ln(two)
    dstar_sign_checks = {
        "log_x_min_gt_log_2": log_x_min[0] > log_two[1],
        "epsilon_below_one": dstar_at_min["epsilon"][1] < 1,
        "q_positive": dstar_at_min["q"][0] > 0,
        "J_lower_positive": dstar_at_min["J_lower"][0] > 0,
        "scaled_prime_square_sum_positive": dstar_at_min[
            "scaled_prime_square_sum_lower"
        ][0]
        > 0,
        "alpha_positive": dstar_at_min["alpha"][0] > 0,
        "epsilon_A_decreases": arithmetic.sub(
            log_two, arithmetic.mul(two, log_x_min)
        )[1]
        < 0,
        "L_over_ell_squared_decreases": arithmetic.sub(log_two, log_x_min)[1]
        < 0,
        "independent_endpoint_overlap": dstar_endpoint_overlap,
    }
    if not all(dstar_sign_checks.values()):
        raise AssertionError(f"D_* monotonicity checks failed: {dstar_sign_checks}")

    boundaries = segment_boundaries()
    segment_rows: list[dict[str, Any]] = []
    minimum_clearance: Interval | None = None
    minimum_segment_index = -1
    maximum_segment_loss: Interval = arithmetic.point(0)
    for index, (left_x, right_x) in enumerate(zip(boundaries, boundaries[1:])):
        left = pointwise_components(
            arithmetic, left_x, pi, zero_sum, log_two_pi, tail_sums
        )
        right_high = sharpened_high_zero_loss(arithmetic, right_x, tail_sums)
        segment_loss = arithmetic.add(
            arithmetic.add(left["low_zero_loss"], left["b2_loss"]), right_high
        )
        segment_clearance = arithmetic.sub(left["dstar"], segment_loss)
        if segment_clearance[0] <= 0:
            raise AssertionError(
                f"nonpositive segment clearance on [{left_x}, {right_x}]"
            )
        if minimum_clearance is None or segment_clearance[0] < minimum_clearance[0]:
            minimum_clearance = segment_clearance
            minimum_segment_index = index
        if segment_loss[1] > maximum_segment_loss[1]:
            maximum_segment_loss = segment_loss
        segment_rows.append(
            {
                "index": index,
                "left": left_x,
                "right": right_x,
                "D_star_lower_at_left": interval_row(left["dstar"]),
                "low_zero_loss_at_left": interval_row(left["low_zero_loss"]),
                "B2_loss_at_left": interval_row(left["b2_loss"]),
                "high_zero_loss_at_right": interval_row(right_high),
                "segment_residual_loss": interval_row(segment_loss),
                "combined_clearance": interval_row(segment_clearance),
            }
        )
    if minimum_clearance is None:
        raise AssertionError("no segments were generated")

    endpoint = pointwise_components(
        arithmetic, X_MAX, pi, zero_sum, log_two_pi, tail_sums
    )
    barrier_low = pointwise_components(
        arithmetic, BARRIER_LOW, pi, zero_sum, log_two_pi, tail_sums
    )
    barrier_high = pointwise_components(
        arithmetic, BARRIER_HIGH, pi, zero_sum, log_two_pi, tail_sums
    )

    right_z = partial_rh_z_data(arithmetic, X_MAX, pi)
    e_interval = arithmetic.exp(one)
    domain_checks = {
        "Dusart_x_min_at_least_3275": X_MIN >= 3_275,
        "Buethe_2016_equation_7_4_x_min_at_least_5000": X_MIN >= 5_000,
        "Buethe_2018_layer_cutoff_domain_x_at_most_1e19": maximum_layer_cutoff[1]
        < Decimal("1e19"),
        "Platt_Trudgian_height_covers_Buethe_requirement": right_z[
            "height_requirement"
        ][1]
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
            "largest layer cutoff", maximum_layer_cutoff, Decimal("1e19")
        ),
        require_upper_below(
            "partial-RH height requirement at X_MAX",
            right_z["height_requirement"],
            Decimal(ZERO_HEIGHT),
        ),
        require_lower_above(
            "sharpened high-zero monotonicity margin",
            high_monotonicity_margin,
            Decimal(0),
        ),
        require_lower_above(
            "delta monotonicity margin", delta_monotonicity_margin, Decimal(0)
        ),
        require_lower_above(
            "B2 monotonicity margin", b2_monotonicity_margin, Decimal(0)
        ),
        require_upper_below(
            "delta at X_MIN is below one", left_z["delta"], Decimal(1)
        ),
        require_lower_above(
            "A_minus at X_MIN is positive", left_z["a_minus"], Decimal(0)
        ),
        require_lower_above(
            "minimum segmented combined clearance",
            minimum_clearance,
            SEGMENT_CLEARANCE_TARGET,
        ),
        require_lower_above(
            "pointwise scalar-envelope clearance at lower barrier",
            barrier_low["clearance"],
            Decimal(0),
        ),
        require_upper_below(
            "pointwise scalar-envelope clearance at upper barrier",
            barrier_high["clearance"],
            Decimal(0),
        ),
    ]

    script_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return {
        "status": "PASS",
        "claim": {
            "range": f"{X_MIN} <= x <= {X_MAX}",
            "profiles": "every CA exponent profile with exact prime support x",
            "combined_theorem": (
                "sqrt(x)log(x)G(n,x)>2e-7 on every segment, using the "
                "increasing analytic D_* lower envelope at the segment left endpoint"
            ),
            "minimum_directed_clearance": interval_row(minimum_clearance),
        },
        "method_improvements_over_1e21_certificate": {
            "two_integrations_by_parts": (
                "For a=1-rho, |J_rho/rho| is bounded by "
                "x^(beta-1)*(h(L)/gamma^2+2*(-h'(L))/gamma^3)."
            ),
            "functional_equation_pairing": (
                "The positive-ordinate multiset is invariant under beta -> 1-beta. "
                "Convexity gives sum x^(beta-1)w_gamma <= "
                "(1+x^-1)/2 sum w_gamma; conjugates then double this bound."
            ),
            "direct_HSW_Stieltjes_tail": (
                "The full HSW main term and error are integrated directly, replacing "
                "the coarser N(t)<=t log(t)/(2pi) tail."
            ),
            "segmented_increasing_D_star": (
                "Each segment uses the increasing closed D_* lower expression at its "
                "left endpoint, rather than the global constant 0.5 or 0.516114."
            ),
        },
        "independent_D_star_certificate": {
            "status": dstar_report["status"],
            "report_path": str(D_STAR_REPORT),
            "report_sha256": dstar_certificate["report_sha256"],
            "certifier_path": str(D_STAR_CERTIFIER),
            "certifier_script_sha256": dstar_certificate["script_sha256"],
            "endpoint_reported": interval_row(reported_dstar_min),
            "endpoint_recomputed": interval_row(
                dstar_at_min["scaled_dstar_lower"]
            ),
            "endpoint_overlap": dstar_endpoint_overlap,
            "variable_q_refinement": (
                "The extended certificate retains q(x)=1-sqrt(2/x), which increases. "
                "The earlier uniform proof froze q at q(X_MIN); both agree at X_MIN."
            ),
            "monotonicity_sign_checks": dstar_sign_checks,
        },
        "tail_zero_sums": {
            "T": ZERO_HEIGHT,
            "S2_positive_ordinates_upper": interval_row(tail_sums["S2_upper"]),
            "S3_positive_ordinates_upper": interval_row(tail_sums["S3_upper"]),
            "HSW_main_contribution": interval_row(tail_sums["main"]),
            "HSW_error_contribution": interval_row(tail_sums["error"]),
        },
        "coverage": {
            "boundary_count": len(boundaries),
            "segment_count": len(segment_rows),
            "coarse_boundary_count": 12,
            "fine_step": X_FINE_STEP,
            "minimum_segment_index": minimum_segment_index,
            "minimum_segment": segment_rows[minimum_segment_index],
            "maximum_segment_residual_loss": interval_row(maximum_segment_loss),
            "segments": segment_rows,
        },
        "endpoint": {
            "x": X_MAX,
            **{name: interval_row(value) for name, value in endpoint.items()},
        },
        "scalar_envelope_barrier_audit": {
            "interpretation": (
                "The same-point scalar lower envelope changes sign in this bracket. "
                "This is a limitation of the present absolute envelope, not a "
                "counterexample and not an impossibility theorem for stronger methods."
            ),
            "positive_point": BARRIER_LOW,
            "positive_clearance": interval_row(barrier_low["clearance"]),
            "negative_point": BARRIER_HIGH,
            "negative_clearance": interval_row(barrier_high["clearance"]),
        },
        "analytic_lemmas": {
            "low_and_B2_endpoint_reduction": (
                "The previous upward certificate proves that the low-zero and B2 "
                "losses decrease for x>=X_MIN."
            ),
            "high_zero_endpoint_reduction": (
                "After scaling, the sharpened high-zero bound is "
                "(exp(L/2)+exp(-L/2))*A(L).  The certified derivative margin is "
                "positive at X_MIN; its positive term increases and its subtracted "
                "rational majorant decreases."
            ),
            "D_star_endpoint_reduction": (
                "epsilon decreases; q, J, (1-epsilon)J, -epsilon*A, -x^-1/2, "
                "and alpha increase.  Thus the displayed D_* lower expression "
                "increases and can be evaluated at each segment left endpoint."
            ),
            "segment_cover": (
                "On [a,b], use D_* at a, low-zero/B2 losses at a, and the "
                "high-zero loss at b.  Every directed segment clearance is positive."
            ),
        },
        "constants": {
            "pi": interval_row(pi),
            "euler_gamma": interval_row(gamma),
            "exact_zero_sum": interval_row(zero_sum),
            "Q_at_X_MIN": interval_row(q_min),
            "Q_at_X_MAX": interval_row(q_max),
            "normalized_layer_majorant": interval_row(normalized_layer_majorant),
            "global_H_coefficient": interval_row(rosser_h_coefficient),
            "maximum_layer_cutoff": interval_row(maximum_layer_cutoff),
            "partial_RH_height_requirement_at_X_MAX": interval_row(
                right_z["height_requirement"]
            ),
            "high_zero_monotonicity_margin": interval_row(
                high_monotonicity_margin
            ),
            "delta_monotonicity_margin": interval_row(delta_monotonicity_margin),
            "B2_monotonicity_margin": interval_row(b2_monotonicity_margin),
        },
        "external_theorem_inputs": [
            "Dusart (1999): next-prime interval and |theta(t)-t|<3.965t/log(t)^2",
            "Rosser-Schoenfeld (1962): theta(t)<1.01624t",
            "Buethe (2016), equation (7.4): partial-RH theta bound",
            "Buethe (2018): theta(u)<u through 1e19, optional layer audit",
            "Platt-Trudgian (2021): RH verified through height 3e12",
            "Hasanalizade-Shen-Wong (2022), Corollary 1.2: explicit N(t) error",
            "Standard zeta functional-equation and conjugation symmetries of zeros",
        ],
        "checks": checks,
        "exact_domain_checks": domain_checks,
        "rigor_scope": {
            "proved_here": (
                "The sharpened high-zero tail, its monotonicity, the variable-q D_* "
                "endpoint expression, all segment comparisons, and domain checks."
            ),
            "trusted": (
                "The published theorem inputs, the earlier hash-validated D_* "
                "certificate, and CPython Decimal transcendental rounding semantics."
            ),
            "remaining_gap": (
                "This covers CA exponent profiles, not the reduction to every integer, "
                "and it is an executable interval certificate rather than a "
                "proof-assistant formalization."
            ),
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "implementation": platform.python_implementation(),
            "decimal_module": decimal.__file__,
            "libmpdec": getattr(decimal, "__libmpdec_version__", None),
            "precision": precision,
            "gamma_harmonic_terms": gamma_terms,
            "script_sha256": script_hash,
            "shared_interval_backend_sha256": hashlib.sha256(
                BASE_CERTIFIER.read_bytes()
            ).hexdigest(),
            "previous_upward_certifier_sha256": hashlib.sha256(
                PREVIOUS_UPWARD_CERTIFIER.read_bytes()
            ).hexdigest(),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Certify the combined CA inequality from 56,048,351 through "
            "1.64967e23 with a sharpened high-zero tail."
        )
    )
    parser.add_argument("--precision", type=int, default=60)
    parser.add_argument("--gamma-terms", type=int, default=100_000)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    certificate = build_certificate(args.precision, args.gamma_terms)
    args.report.write_text(
        json.dumps(certificate, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    minimum = certificate["claim"]["minimum_directed_clearance"]
    barrier = certificate["scalar_envelope_barrier_audit"]
    print(f"status: {certificate['status']}")
    print(f"range: {certificate['claim']['range']}")
    print(f"segment count: {certificate['coverage']['segment_count']}")
    print(f"minimum combined clearance lower: {minimum['lower']}")
    print(
        "partial-RH requirement upper: "
        f"{certificate['constants']['partial_RH_height_requirement_at_X_MAX']['upper']}"
    )
    print(
        "scalar-envelope barrier: "
        f"positive at {barrier['positive_point']}, negative at {barrier['negative_point']}"
    )
    print(f"wrote {args.report}")


if __name__ == "__main__":
    main()
