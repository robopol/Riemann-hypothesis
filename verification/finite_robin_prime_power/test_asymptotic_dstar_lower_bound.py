from __future__ import annotations

import decimal
import hashlib
import json
import platform
import sys
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_REPORT = ROOT / "asymptotic_dstar_lower_bound_report.json"

X0 = 56_048_351
ETA = Decimal("3.965")
TARGET = Decimal("0.5")

Interval = tuple[Decimal, Decimal]


class DecimalIntervals:
    """Small outward-rounded interval layer for the endpoint audit."""

    def __init__(self, precision: int = 60) -> None:
        self.precision = precision
        self.nearest = Context(prec=precision, rounding=ROUND_HALF_EVEN)
        self.floor = Context(prec=precision, rounding=ROUND_FLOOR)
        self.ceiling = Context(prec=precision, rounding=ROUND_CEILING)
        self.padding = Context(prec=precision + 8, rounding=ROUND_HALF_EVEN)

    def point(self, value: int | str | Decimal) -> Interval:
        number = value if isinstance(value, Decimal) else Decimal(value)
        return number, number

    def _ulp(self, value: Decimal) -> Decimal:
        return Decimal(1).scaleb(value.adjusted() - self.precision + 1)

    def _padded_unary(self, operation: str, value: Decimal) -> Interval:
        if operation == "ln":
            midpoint = self.nearest.ln(value)
        elif operation == "sqrt":
            midpoint = self.nearest.sqrt(value)
        else:
            raise ValueError(f"unsupported unary operation: {operation}")
        ulp = self._ulp(midpoint)
        return (
            self.padding.subtract(midpoint, ulp),
            self.padding.add(midpoint, ulp),
        )

    def ln(self, value: Interval) -> Interval:
        if value[0] <= 0:
            raise ValueError("logarithm interval must be positive")
        return (
            self._padded_unary("ln", value[0])[0],
            self._padded_unary("ln", value[1])[1],
        )

    def sqrt(self, value: Interval) -> Interval:
        if value[0] < 0:
            raise ValueError("square-root interval must be nonnegative")
        return (
            self._padded_unary("sqrt", value[0])[0],
            self._padded_unary("sqrt", value[1])[1],
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
        low = [self.floor.multiply(a, b) for a in left for b in right]
        high = [self.ceiling.multiply(a, b) for a in left for b in right]
        return min(low), max(high)

    def div(self, numerator: Interval, denominator: Interval) -> Interval:
        if denominator[0] <= 0 <= denominator[1]:
            raise ZeroDivisionError("denominator interval contains zero")
        low = [self.floor.divide(a, b) for a in numerator for b in denominator]
        high = [self.ceiling.divide(a, b) for a in numerator for b in denominator]
        return min(low), max(high)


def interval_row(value: Interval) -> dict[str, str]:
    return {"lower": str(value[0]), "upper": str(value[1])}


def compute_endpoint_certificate(precision: int = 60) -> dict[str, object]:
    arithmetic = DecimalIntervals(precision)
    one = arithmetic.point(1)
    two = arithmetic.point(2)
    three = arithmetic.point(3)
    half = arithmetic.div(one, two)
    x = arithmetic.point(X0)
    log_x = arithmetic.ln(x)
    log_two = arithmetic.ln(two)
    ell = arithmetic.div(arithmetic.add(log_x, log_two), two)
    sqrt_two = arithmetic.sqrt(two)
    sqrt_x = arithmetic.sqrt(x)
    y = arithmetic.sqrt(arithmetic.mul(two, x))
    r = arithmetic.sqrt(arithmetic.div(x, two))

    ell_squared = arithmetic.mul(ell, ell)
    epsilon = arithmetic.div(arithmetic.point(ETA), ell_squared)
    one_minus_epsilon = arithmetic.sub(one, epsilon)

    a_term = arithmetic.div(
        log_x,
        arithmetic.mul(sqrt_two, ell),
    )
    q0 = arithmetic.sub(one, arithmetic.div(one, r))
    j_inner = arithmetic.sub(
        arithmetic.div(q0, ell),
        arithmetic.div(one, ell_squared),
    )
    j_lower = arithmetic.mul(
        arithmetic.div(log_x, sqrt_two),
        j_inner,
    )

    scaled_s2_lower = arithmetic.sub(
        arithmetic.sub(
            arithmetic.mul(one_minus_epsilon, j_lower),
            arithmetic.mul(
                two,
                arithmetic.mul(epsilon, a_term),
            ),
        ),
        arithmetic.div(one, sqrt_x),
    )
    alpha = arithmetic.sub(
        half,
        arithmetic.div(one, arithmetic.mul(three, y)),
    )
    scaled_dstar_lower = arithmetic.mul(alpha, scaled_s2_lower)
    margin = arithmetic.sub(scaled_dstar_lower, arithmetic.point(TARGET))

    monotonicity_checks = {
        "log_x0_gt_log_2": log_x[0] > log_two[1],
        "epsilon_below_one": epsilon[1] < 1,
        "j_lower_positive": j_lower[0] > 0,
        "scaled_s2_lower_positive": scaled_s2_lower[0] > 0,
        "epsilon_times_A_derivative_numerator_negative": (
            arithmetic.sub(log_two, arithmetic.mul(two, log_x))[1] < 0
        ),
        "J_derivative_positive_terms": (
            q0[0] > 0 and arithmetic.sub(log_x, log_two)[0] > 0
        ),
        "alpha_positive": alpha[0] > 0,
        "strict_target_pass": scaled_dstar_lower[0] > TARGET,
    }

    return {
        "claim": (
            "For every real x >= 56,048,351, the exact full-support envelope "
            "satisfies sqrt(x) log(x) D_*(x) > 0.5."
        ),
        "status": "PASS" if all(monotonicity_checks.values()) else "FAIL",
        "method": {
            "prime_subset": "sqrt(2x) < p <= x",
            "per_prime_lower_bound": (
                "D_p=1/p-log(1+1/p) > "
                "(1/2-1/(3 sqrt(2x)))/p^2"
            ),
            "theta_input": (
                "Dusart 1999, Theorem 4: |theta(t)-t| < "
                "3.965 t/log(t)^2 for t>=2"
            ),
            "uniformity": (
                "The closed lower expression is increasing in log(x); "
                "the report evaluates its left endpoint."
            ),
        },
        "endpoint": {
            "x0": X0,
            "log_x": interval_row(log_x),
            "ell_log_sqrt_2x": interval_row(ell),
            "y_sqrt_2x": interval_row(y),
            "r_sqrt_x_over_2": interval_row(r),
            "epsilon": interval_row(epsilon),
            "A": interval_row(a_term),
            "q0": interval_row(q0),
            "J_lower": interval_row(j_lower),
            "scaled_prime_square_sum_lower": interval_row(scaled_s2_lower),
            "per_prime_factor_alpha": interval_row(alpha),
            "scaled_dstar_lower": interval_row(scaled_dstar_lower),
            "margin_over_target": interval_row(margin),
        },
        "monotonicity_checks": monotonicity_checks,
        "parameters": {
            "precision": precision,
            "dusart_eta": str(ETA),
            "target": str(TARGET),
        },
    }


def main() -> None:
    report = compute_endpoint_certificate()
    script_path = Path(__file__).resolve()
    report["runtime_environment"] = {
        "python": sys.version,
        "platform": platform.platform(),
        "decimal_module": decimal.__file__,
        "libmpdec_version": getattr(decimal, "__libmpdec_version__", "unknown"),
        "certifier_script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    DEFAULT_REPORT.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    endpoint = report["endpoint"]
    print(f"wrote {DEFAULT_REPORT}")
    print(
        f"status={report['status']} "
        f"lower={endpoint['scaled_dstar_lower']['lower']} "
        f"margin={endpoint['margin_over_target']['lower']}"
    )


if __name__ == "__main__":
    main()
