from __future__ import annotations

import argparse
import decimal
import hashlib
import json
import platform
import sys
from typing import Any
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_REPORT = ROOT / "ca_residual_finite_certificate_report.json"
DEFAULT_BUFFER_REPORT = ROOT / "ca_exact_buffer_interval_certificate_report.json"
DEFAULT_BUFFER_SCRIPT = ROOT / "test_ca_exact_buffer_interval_certifier.py"

X_MIN = 3_329_267
X_MAX = 56_048_351
ZERO_HEIGHT = 3_000_000_000_000
BUFFER_TARGET = Decimal("0.7825")
RESIDUAL_LOSS_TARGET = Decimal("0.0585")

Interval = tuple[Decimal, Decimal]


class DecimalIntervals:
    """Small outward-rounded Decimal interval backend for scalar constants."""

    def __init__(self, precision: int) -> None:
        if precision < 40:
            raise ValueError("precision must be at least 40 decimal digits")
        self.precision = precision
        self.nearest = Context(prec=precision, rounding=ROUND_HALF_EVEN)
        self.floor = Context(prec=precision, rounding=ROUND_FLOOR)
        self.ceiling = Context(prec=precision, rounding=ROUND_CEILING)
        self.padding = Context(prec=precision + 10, rounding=ROUND_HALF_EVEN)

    def point(self, value: int | str | Decimal) -> Interval:
        number = value if isinstance(value, Decimal) else Decimal(value)
        return number, number

    def _ulp(self, value: Decimal) -> Decimal:
        return Decimal(1).scaleb(value.adjusted() - self.precision + 1)

    def _padded_unary(self, operation: str, value: Decimal) -> Interval:
        if operation == "ln":
            midpoint = self.nearest.ln(value)
        elif operation == "exp":
            midpoint = self.nearest.exp(value)
        elif operation == "sqrt":
            midpoint = self.nearest.sqrt(value)
        else:
            raise ValueError(f"unsupported unary operation {operation!r}")
        padding = self.padding.multiply(Decimal(2), self._ulp(midpoint))
        return (
            self.padding.subtract(midpoint, padding),
            self.padding.add(midpoint, padding),
        )

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
        lower = [self.floor.multiply(a, b) for a in left for b in right]
        upper = [self.ceiling.multiply(a, b) for a in left for b in right]
        return min(lower), max(upper)

    def div(self, numerator: Interval, denominator: Interval) -> Interval:
        if denominator[0] <= 0 <= denominator[1]:
            raise ZeroDivisionError("division interval contains zero")
        lower = [self.floor.divide(a, b) for a in numerator for b in denominator]
        upper = [self.ceiling.divide(a, b) for a in numerator for b in denominator]
        return min(lower), max(upper)

    def ln(self, value: Interval) -> Interval:
        if value[0] <= 0:
            raise ValueError("logarithm interval must be positive")
        return (
            self._padded_unary("ln", value[0])[0],
            self._padded_unary("ln", value[1])[1],
        )

    def exp(self, value: Interval) -> Interval:
        return (
            self._padded_unary("exp", value[0])[0],
            self._padded_unary("exp", value[1])[1],
        )

    def sqrt(self, value: Interval) -> Interval:
        if value[0] < 0:
            raise ValueError("square-root interval must be nonnegative")
        return (
            self._padded_unary("sqrt", value[0])[0],
            self._padded_unary("sqrt", value[1])[1],
        )

    def rational(self, numerator: int, denominator: int) -> Interval:
        return self.div(self.point(numerator), self.point(denominator))

    def rational_power(self, value: Interval, numerator: int, denominator: int) -> Interval:
        exponent = self.rational(numerator, denominator)
        return self.exp(self.mul(exponent, self.ln(value)))


def interval_row(value: Interval) -> dict[str, str]:
    return {"lower": str(value[0]), "upper": str(value[1]), "width": str(value[1] - value[0])}


def alternating_arctangent_reciprocal(
    arithmetic: DecimalIntervals, denominator: int, last_even_index: int = 100
) -> Interval:
    """Enclose atan(1/q) between consecutive alternating partial sums."""
    if last_even_index % 2:
        raise ValueError("last_even_index must be even")
    lower = arithmetic.point(0)
    upper = arithmetic.point(0)
    partial = arithmetic.point(0)
    q = denominator
    for index in range(last_even_index + 1):
        term = arithmetic.rational(1, (2 * index + 1) * (q ** (2 * index + 1)))
        partial = arithmetic.add(partial, term) if index % 2 == 0 else arithmetic.sub(partial, term)
        if index == last_even_index - 1:
            lower = partial
        elif index == last_even_index:
            upper = partial
    return lower[0], upper[1]


def machin_pi_interval(arithmetic: DecimalIntervals) -> Interval:
    atan_five = alternating_arctangent_reciprocal(arithmetic, 5)
    atan_239 = alternating_arctangent_reciprocal(arithmetic, 239)
    return arithmetic.sub(
        arithmetic.mul(arithmetic.point(16), atan_five),
        arithmetic.mul(arithmetic.point(4), atan_239),
    )


def euler_gamma_interval(arithmetic: DecimalIntervals, terms: int) -> Interval:
    """Use 1/(2n+1) < H_n-log(n)-gamma < 1/(2n)."""
    harmonic = arithmetic.point(0)
    for value in range(1, terms + 1):
        harmonic = arithmetic.add(harmonic, arithmetic.rational(1, value))
    log_n = arithmetic.ln(arithmetic.point(terms))
    center = arithmetic.sub(harmonic, log_n)
    lower = arithmetic.sub(center, arithmetic.rational(1, 2 * terms))[0]
    upper = arithmetic.sub(center, arithmetic.rational(1, 2 * terms + 1))[1]
    return lower, upper


def require_upper_below(name: str, value: Interval, target: Decimal) -> dict[str, str]:
    if value[1] >= target:
        raise AssertionError(f"{name}: upper {value[1]} is not below {target}")
    return {
        "name": name,
        "status": "PASS",
        "upper": str(value[1]),
        "target": str(target),
        "margin": str(target - value[1]),
    }


def require_lower_above(name: str, value: Interval, target: Decimal) -> dict[str, str]:
    if value[0] <= target:
        raise AssertionError(f"{name}: lower {value[0]} is not above {target}")
    return {
        "name": name,
        "status": "PASS",
        "lower": str(value[0]),
        "target": str(target),
        "margin": str(value[0] - target),
    }


def load_buffer_certificate(report_path: Path, script_path: Path) -> dict[str, Any]:
    """Validate and summarize the independent full-support D_* certificate."""
    report_bytes = report_path.read_bytes()
    report = json.loads(report_bytes)
    script_hash = hashlib.sha256(script_path.read_bytes()).hexdigest()
    embedded_script_hash = report["runtime_environment"]["certifier_script_sha256"]

    if report["status"] != "PASS":
        raise AssertionError("the independent D_* certificate is not PASS")
    if embedded_script_hash != script_hash:
        raise AssertionError(
            "the D_* report does not match the current certifier script: "
            f"report has {embedded_script_hash}, current script has {script_hash}"
        )

    parameters = report["parameters"]
    counts = report["status_counts"]
    coverage = report["coverage"]
    if parameters["min_prime"] > X_MIN or parameters["max_prime"] < X_MAX:
        raise AssertionError("the D_* certificate does not cover the residual range")
    if Decimal(parameters["target"]) != BUFFER_TARGET:
        raise AssertionError("the D_* certificate used a different strict target")
    if counts["fail"] != 0 or counts["inconclusive"] != 0:
        raise AssertionError("the D_* certificate has failed or inconclusive supports")
    if counts["pass"] != counts["support_count"]:
        raise AssertionError("the D_* certificate did not pass every swept support")
    if coverage["first_support"] != parameters["min_prime"]:
        raise AssertionError("the D_* coverage starts at an unexpected support")
    if coverage["last_support"] != parameters["max_prime"]:
        raise AssertionError("the D_* coverage ends at an unexpected support")

    minimum = report["minimum"]
    minimum_interval = (
        Decimal(minimum["scaled_buffer"]["lower"]),
        Decimal(minimum["scaled_buffer"]["upper"]),
    )
    if minimum_interval[0] <= BUFFER_TARGET:
        raise AssertionError("the certified D_* minimum is not strictly above 0.7825")

    return {
        "report": report,
        "minimum_interval": minimum_interval,
        "report_sha256": hashlib.sha256(report_bytes).hexdigest(),
        "script_sha256": script_hash,
    }


def build_certificate(
    precision: int,
    gamma_terms: int,
    buffer_report_path: Path = DEFAULT_BUFFER_REPORT,
    buffer_script_path: Path = DEFAULT_BUFFER_SCRIPT,
) -> dict[str, object]:
    arithmetic = DecimalIntervals(precision)
    buffer_certificate = load_buffer_certificate(buffer_report_path, buffer_script_path)
    buffer_report = buffer_certificate["report"]
    buffer_minimum = buffer_certificate["minimum_interval"]
    one = arithmetic.point(1)
    two = arithmetic.point(2)
    three = arithmetic.point(3)
    four = arithmetic.point(4)
    ten = arithmetic.point(10)

    pi = machin_pi_interval(arithmetic)
    gamma = euler_gamma_interval(arithmetic, gamma_terms)
    log_two = arithmetic.ln(two)
    log_x_min = arithmetic.ln(arithmetic.point(X_MIN))
    log_x_max = arithmetic.ln(arithmetic.point(X_MAX))
    sqrt_x_min = arithmetic.sqrt(arithmetic.point(X_MIN))
    sqrt_x_max = arithmetic.sqrt(arithmetic.point(X_MAX))

    # Dusart's next-prime theorem gives q+1 <= Q(x).
    log_x_min_squared = arithmetic.mul(log_x_min, log_x_min)
    q_ratio_extra = arithmetic.div(one, arithmetic.mul(two, log_x_min_squared))
    q_upper_min = arithmetic.add(
        arithmetic.mul(arithmetic.point(X_MIN), arithmetic.add(one, q_ratio_extra)),
        one,
    )
    q_over_x_min = arithmetic.div(q_upper_min, arithmetic.point(X_MIN))

    log_x_max_squared = arithmetic.mul(log_x_max, log_x_max)
    q_ratio_extra_max = arithmetic.div(one, arithmetic.mul(two, log_x_max_squared))
    q_upper_max = arithmetic.add(
        arithmetic.mul(arithmetic.point(X_MAX), arithmetic.add(one, q_ratio_extra_max)),
        one,
    )
    maximum_layer_cutoff = arithmetic.sqrt(arithmetic.mul(two, q_upper_max))

    # Analytic majorant for F(Q)=sum_{s>=2,(sQ)^(1/s)>2}(sQ)^(1/s).
    log2_q = arithmetic.div(arithmetic.ln(q_upper_min), log_two)
    twice_log2_q = arithmetic.mul(two, log2_q)
    b_q = arithmetic.add(log2_q, arithmetic.div(arithmetic.ln(twice_log2_q), log_two))
    b_minus_four = arithmetic.sub(b_q, four)
    sqrt_two = arithmetic.sqrt(two)
    term_three = arithmetic.mul(
        arithmetic.rational_power(three, 1, 3),
        arithmetic.rational_power(q_upper_min, -1, 6),
    )
    term_four = arithmetic.mul(
        arithmetic.rational_power(four, 1, 4),
        arithmetic.rational_power(q_upper_min, -1, 4),
    )
    term_five = arithmetic.mul(
        arithmetic.mul(arithmetic.rational_power(arithmetic.point(5), 1, 5), b_minus_four),
        arithmetic.rational_power(q_upper_min, -3, 10),
    )
    normalized_layer_sum_majorant = arithmetic.add(
        arithmetic.add(sqrt_two, term_three), arithmetic.add(term_four, term_five)
    )

    # The last term decreases once 0.3(B-4)>B', with y=log Q.
    log_q = arithmetic.ln(q_upper_min)
    b_derivative = arithmetic.div(
        arithmetic.add(one, arithmetic.div(one, log_q)), log_two
    )
    tail_derivative_margin = arithmetic.sub(
        arithmetic.mul(arithmetic.rational(3, 10), b_minus_four), b_derivative
    )

    two_sqrt_q_over_x = arithmetic.mul(two, arithmetic.sqrt(q_over_x_min))
    global_theta_layer_coefficient = arithmetic.mul(
        arithmetic.point("1.01624"), two_sqrt_q_over_x
    )

    # The finite Büthe theta theorem gives theta(u)<u for all layer cutoffs.
    finite_layer_coefficient = two_sqrt_q_over_x

    # Büthe gives -1.95 sqrt(x) <= theta(x)-x < -0.05 sqrt(x).
    # With 0<=H<2.003 sqrt(x), this implies |z|<2/sqrt(x).
    delta_min = arithmetic.div(two, sqrt_x_min)
    a_min = arithmetic.add(log_x_min, arithmetic.ln(arithmetic.sub(one, delta_min)))
    b2_curvature = arithmetic.div(
        arithmetic.add(a_min, one),
        arithmetic.mul(
            arithmetic.mul(arithmetic.sub(one, delta_min), arithmetic.sub(one, delta_min)),
            arithmetic.mul(a_min, a_min),
        ),
    )
    b2_scaled_loss = arithmetic.mul(
        arithmetic.mul(arithmetic.div(two, sqrt_x_min), log_x_min), b2_curvature
    )

    # Exact zero sum S0=2+gamma-log(4*pi).
    four_pi = arithmetic.mul(four, pi)
    zero_sum = arithmetic.sub(arithmetic.add(two, gamma), arithmetic.ln(four_pi))
    low_kernel_factor = arithmetic.add(
        one,
        arithmetic.add(
            arithmetic.div(three, log_x_min),
            arithmetic.div(four, arithmetic.mul(log_x_min, log_x_min)),
        ),
    )
    log_two_pi = arithmetic.ln(arithmetic.mul(two, pi))
    low_zero_scaled_loss = arithmetic.add(
        arithmetic.mul(zero_sum, low_kernel_factor),
        arithmetic.div(log_two_pi, sqrt_x_min),
    )

    # High zeros: HSW gives N(t)<=t log(t)/(2*pi) for t>=T.
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
    high_zero_scaled_loss = arithmetic.div(
        arithmetic.mul(
            arithmetic.mul(
                arithmetic.mul(four, sqrt_x_max),
                arithmetic.div(arithmetic.add(log_x_max, one), log_x_max),
            ),
            arithmetic.add(log_zero_height, one),
        ),
        arithmetic.mul(pi, zero_height),
    )

    partial_rh_lhs = arithmetic.mul(
        arithmetic.point("4.92"),
        arithmetic.sqrt(arithmetic.div(arithmetic.point(X_MAX), log_x_max)),
    )

    total_residual_loss = arithmetic.add(
        arithmetic.add(low_zero_scaled_loss, high_zero_scaled_loss), b2_scaled_loss
    )
    target_clearance = arithmetic.sub(arithmetic.point(BUFFER_TARGET), total_residual_loss)
    certified_clearance = arithmetic.sub(buffer_minimum, total_residual_loss)

    domain_checks = {
        "Dusart_x_min_at_least_3275": X_MIN >= 3_275,
        "Buethe_theta_lower_domain_x_min_at_least_1423": X_MIN >= 1_423,
        "Buethe_theta_upper_domain_x_max_at_most_1e19": X_MAX <= 10**19,
        "Buethe_partial_RH_equation_7_4_x_min_at_least_5000": X_MIN >= 5_000,
        "HSW_zero_count_height_at_least_e": Decimal(ZERO_HEIGHT) > e_interval[1],
    }
    failed_domains = [name for name, passed in domain_checks.items() if not passed]
    if failed_domains:
        raise AssertionError(f"external theorem domain checks failed: {failed_domains}")

    checks = [
        require_upper_below(
            "normalized layer-sum majorant at Q(X_MIN)",
            normalized_layer_sum_majorant,
            Decimal(2),
        ),
        require_lower_above(
            "tail-majorant monotonicity margin",
            tail_derivative_margin,
            Decimal(0),
        ),
        require_upper_below(
            "global H/sqrt(x) coefficient",
            global_theta_layer_coefficient,
            Decimal("2.035"),
        ),
        require_upper_below(
            "finite H/sqrt(x) coefficient",
            finite_layer_coefficient,
            Decimal("2.003"),
        ),
        require_upper_below(
            "largest possible layer cutoff is inside Buethe finite domain",
            maximum_layer_cutoff,
            Decimal("1e19"),
        ),
        require_upper_below(
            "scaled B2 loss",
            b2_scaled_loss,
            Decimal("0.00118"),
        ),
        require_lower_above(
            "zero-count upper-bound margin at T",
            zero_count_margin_at_t,
            Decimal(0),
        ),
        require_lower_above(
            "zero-count margin derivative for t>=T",
            zero_count_derivative_margin,
            Decimal(0),
        ),
        require_upper_below(
            "high-zero scaled loss",
            high_zero_scaled_loss,
            Decimal("0.0000001"),
        ),
        require_upper_below(
            "Buethe partial-RH height requirement",
            partial_rh_lhs,
            Decimal(ZERO_HEIGHT),
        ),
        require_upper_below(
            "total scaled residual loss",
            total_residual_loss,
            RESIDUAL_LOSS_TARGET,
        ),
        require_lower_above(
            "scaled Robin clearance from the rounded D_* target",
            target_clearance,
            Decimal("0.724083229312263"),
        ),
        require_lower_above(
            "scaled Robin clearance from the certified D_* minimum",
            certified_clearance,
            Decimal("0.72418"),
        ),
    ]

    script_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return {
        "status": "PASS",
        "claim": {
            "range": f"every CA profile with prime support {X_MIN} <= x <= {X_MAX}",
            "residual_bound": (
                "sqrt(x)log(x)*(R_infinity(x)+B2(n,x)) > -0.0585"
            ),
            "combined_finite_CA_theorem": (
                "the independent full-support interval certificate proves "
                "sqrt(x)log(x)D_*(x)>0.7825, hence "
                "sqrt(x)log(x)G(n,x)>0.724083229312263"
            ),
        },
        "independent_D_star_certificate": {
            "status": buffer_report["status"],
            "report_path": str(buffer_report_path),
            "report_sha256": buffer_certificate["report_sha256"],
            "certifier_path": str(buffer_script_path),
            "certifier_script_sha256": buffer_certificate["script_sha256"],
            "embedded_certifier_script_sha256": buffer_report["runtime_environment"][
                "certifier_script_sha256"
            ],
            "prime_stream_sha256_little_endian_uint64": buffer_report["coverage"][
                "prime_stream_sha256_little_endian_uint64"
            ],
            "swept_range": {
                "first_support": buffer_report["coverage"]["first_support"],
                "last_support": buffer_report["coverage"]["last_support"],
                "support_count": buffer_report["status_counts"]["support_count"],
            },
            "minimum_support": buffer_report["minimum"]["support_prime"],
            "minimum_scaled_D_star": interval_row(buffer_minimum),
            "strict_target": str(BUFFER_TARGET),
        },
        "external_theorem_inputs": [
            {
                "source": "Dusart (1999), prime in a short interval",
                "input": (
                    "for x>=3275, the next prime q satisfies "
                    "q<=x(1+1/(2 log(x)^2))"
                ),
                "used_for": "uniform tie-safe CA layer cutoff",
            },
            {
                "source": "Rosser-Schoenfeld (1962), equation (3.24)",
                "input": "theta(t)<1.01624 t for every t>0",
                "used_for": "standalone global H<2.035 sqrt(x) lemma",
            },
            {
                "source": "Buethe, Math. Comp. 87 (2018), Theorem 2",
                "input": (
                    "x-theta(x)<=1.95 sqrt(x) for 1423<=x<=1e19 and "
                    "x-theta(x)>0.05 sqrt(x) for 1<=x<=1e19"
                ),
                "used_for": "finite theta sign, H<2.003 sqrt(x), and |z|<2/sqrt(x)",
            },
            {
                "source": "Platt-Trudgian, Bull. LMS 53 (2021)",
                "input": "all zeta zeros with 0<Im(rho)<=3e12 lie on Re(rho)=1/2",
                "used_for": "low-zero kernel bound",
            },
            {
                "source": "Hasanalizade-Shen-Wong, JNT 235 (2022), Corollary 1.2",
                "input": (
                    "|N(t)-t/(2pi)log(t/(2pi e))| <= "
                    "0.1038log(t)+0.2573loglog(t)+9.3675 for t>=e"
                ),
                "used_for": "unverified high-zero tail",
            },
            {
                "source": "Buethe, Math. Comp. 85 (2016), Theorem 2 and equation (7.4)",
                "input": (
                    "partial RH through T and 4.92sqrt(x/log x)<=T imply the "
                    "Schoenfeld theta bound; equation (7.4) gives the stronger "
                    "sqrt(x)log(x)(log(x)-2)/(8pi) form for x>=5000"
                ),
                "used_for": "audited fallback; not needed for the final sharper B2 constant",
            },
        ],
        "analytic_lemmas": {
            "tie_safe_layer": (
                "If a_p>=s>=2, then tau_1(q)<=epsilon<=tau_s(p), "
                "so p^s log p<(q+1)log q and p<(s(q+1))^(1/s)."
            ),
            "layer_sum": (
                "F(Q)/sqrt(Q) is bounded by sqrt(2)+3^(1/3)Q^(-1/6)+"
                "4^(1/4)Q^(-1/4)+5^(1/5)(B(Q)-4)Q^(-3/10), "
                "B(Q)=log_2(Q)+log_2(2log_2(Q)); this majorant decreases."
            ),
            "b2_taylor": (
                "b_L''(z)=-(A+1)/((1+z)^2 A^2), A=L+log(1+z); "
                "Taylor at z=0 and |z|<2/sqrt(x) give the stated loss."
            ),
            "low_zero_kernel": (
                "integration by parts gives |J_rho|<=x^-1/2/|1-rho|*"
                "[(L+1)/L^2+2(L+2)/L^3]; the exact zero sum is "
                "2+gamma-log(4pi)."
            ),
            "high_zero_tail": (
                "N(t)<=tlog(t)/(2pi) for t>=T implies "
                "sum_{gamma>T}gamma^-2<=(log(T)+1)/(pi T)."
            ),
            "termwise_zero_integration": (
                "Integrate symmetric finite zero truncations first. After integration, "
                "|J_rho/rho|<=2h(L)/gamma^2 for high zeros, and N(t)=O(t log t) "
                "makes the integrated zero series absolutely convergent; compact "
                "cutoffs and their limits therefore justify the interchange."
            ),
            "monotonic_extrema": [
                "Q(x)/x=1+1/(2log(x)^2)+1/x decreases",
                "the normalized low-zero loss decreases with x",
                "the high-zero loss increases on the finite interval",
                (
                    "the B2 majorant is bounded by an explicitly decreasing "
                    "e^(-L/2) times a decreasing rational function of L"
                ),
            ],
        },
        "constants": {
            "pi": interval_row(pi),
            "euler_gamma": interval_row(gamma),
            "exact_zero_sum": interval_row(zero_sum),
            "Q_at_x_min": interval_row(q_upper_min),
            "Q_at_x_max": interval_row(q_upper_max),
            "maximum_layer_cutoff": interval_row(maximum_layer_cutoff),
            "normalized_layer_sum_majorant": interval_row(normalized_layer_sum_majorant),
            "global_H_coefficient": interval_row(global_theta_layer_coefficient),
            "finite_H_coefficient": interval_row(finite_layer_coefficient),
            "delta_at_x_min": interval_row(delta_min),
            "scaled_B2_loss": interval_row(b2_scaled_loss),
            "scaled_low_zero_loss": interval_row(low_zero_scaled_loss),
            "scaled_high_zero_loss": interval_row(high_zero_scaled_loss),
            "scaled_total_residual_loss": interval_row(total_residual_loss),
            "scaled_clearance_from_rounded_D_star_target": interval_row(target_clearance),
            "scaled_clearance_from_certified_D_star_minimum": interval_row(
                certified_clearance
            ),
            "partial_RH_requirement_at_x_max": interval_row(partial_rh_lhs),
            "zero_count_margin_at_T": interval_row(zero_count_margin_at_t),
            "zero_count_derivative_margin": interval_row(zero_count_derivative_margin),
        },
        "checks": checks,
        "exact_domain_checks": domain_checks,
        "rigor_scope": {
            "proved_here": (
                "All algebraic reductions, monotonic endpoint reductions, and scalar "
                "constant comparisons listed in analytic_lemmas."
            ),
            "external": (
                "The published theorems listed in external_theorem_inputs are accepted "
                "as inputs and are not recomputed by this script."
            ),
            "combined_status": (
                "The separate interval certificate passed every swept prime support, "
                "so the combined finite CA statement is complete subject to the "
                "published theorem inputs and executable-certificate trust assumptions."
            ),
            "remaining_gap": (
                "This is not a proof-assistant formalization. The external published "
                "theorems, both Decimal runtimes, deterministic sieve, and elementary "
                "certificate lemmas remain trusted inputs. Reduction from CA supports "
                "to every integer and any asymptotic extension are separate problems."
            ),
            "decimal_assumption": (
                "CPython Decimal ln, exp, and sqrt have documented correct-rounding "
                "semantics; two ulps are added to every transcendental endpoint."
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
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Certify scalar constants in the finite CA residual lemma."
    )
    parser.add_argument("--precision", type=int, default=60)
    parser.add_argument("--gamma-terms", type=int, default=100_000)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--buffer-report", type=Path, default=DEFAULT_BUFFER_REPORT)
    parser.add_argument("--buffer-script", type=Path, default=DEFAULT_BUFFER_SCRIPT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    certificate = build_certificate(
        args.precision,
        args.gamma_terms,
        args.buffer_report,
        args.buffer_script,
    )
    args.report.write_text(json.dumps(certificate, indent=2), encoding="utf-8")
    constants = certificate["constants"]
    print(f"status: {certificate['status']}")
    print(f"scaled B2 loss upper: {constants['scaled_B2_loss']['upper']}")
    print(f"scaled low-zero loss upper: {constants['scaled_low_zero_loss']['upper']}")
    print(f"scaled high-zero loss upper: {constants['scaled_high_zero_loss']['upper']}")
    print(f"scaled total loss upper: {constants['scaled_total_residual_loss']['upper']}")
    print(
        "scaled clearance lower from rounded target: "
        f"{constants['scaled_clearance_from_rounded_D_star_target']['lower']}"
    )
    print(
        "scaled clearance lower from certified D_* minimum: "
        f"{constants['scaled_clearance_from_certified_D_star_minimum']['lower']}"
    )
    print(f"wrote {args.report}")


if __name__ == "__main__":
    main()
