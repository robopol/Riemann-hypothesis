#!/usr/bin/env python3
"""Regression checks for the exploratory prime-power B-cell dynamics.

These checks validate algebraic/numerical consistency on a moderate finite
range. They are not an interval certificate, a polynomial-envelope proof, or
a proof of the Riemann hypothesis.
"""

from __future__ import annotations

import argparse
import math

from exploratory_b_cell_dynamics import (
    build_report,
    cell_b_and_derivative,
    forcing_g_and_derivative,
)


def check_open_cell_ode() -> None:
    """Check the smooth ODE independently by differentiating B' numerically."""

    t = 1.1
    s_lambda = math.log(2.0) / 2.0
    psi = math.log(2.0)
    step = 1.0e-5
    b_value = float(cell_b_and_derivative(t, s_lambda, psi)[0])
    derivative_minus = float(
        cell_b_and_derivative(t - step, s_lambda, psi)[1]
    )
    derivative_plus = float(
        cell_b_and_derivative(t + step, s_lambda, psi)[1]
    )
    second_derivative = (derivative_plus - derivative_minus) / (2.0 * step)
    forcing = float(forcing_g_and_derivative(t)[0])
    residual = second_derivative - 0.25 * b_value + forcing
    if abs(residual) >= 2.0e-9:
        raise AssertionError(f"open-cell ODE residual is too large: {residual}")


def check_event_scan() -> None:
    """Run a finite event scan and check every implemented exact identity."""

    report = build_report(
        argparse.Namespace(
            limit=200_000,
            scan_min=2,
            bin_width=1.0,
            powers=[0.0, 1.0, 2.0],
        )
    )
    if not report["checks"]["all_pass"]:
        raise AssertionError(f"report consistency checks failed: {report['checks']}")

    directions = report["stationary_extrema"][
        "stationary_derivative_direction_counts"
    ]
    if directions["negative_to_positive"] != 0:
        raise AssertionError(f"unexpected interior minimum: {directions}")

    cone = report["Phi_first_crossing_energy"]["negative_H_cone_margins"]
    if not cone["both_one_sided_states_strictly_inside_cone"]:
        raise AssertionError("the finite first-crossing cone check failed")
    if cone["minimum_pre_kick_ratio_R_plus_Phi_prime_left_over_j"] <= 1.0:
        raise AssertionError(f"the finite kick inequality failed: {cone}")

    feasibility = report["constant_quadratic_jump_feasibility"]
    if feasibility[
        "constant_positive_definite_quadratic_can_be_nonincreasing_at_all_scanned_jumps"
    ]:
        raise AssertionError("the disqualified constant quadratic family reappeared")

    energy_rows = report["regularized_Bhat_polynomial_energy"]["rows"]
    if max(row["max_abs_exact_jump_formula_error"] for row in energy_rows) >= 1e-8:
        raise AssertionError("a regularized polynomial-energy jump identity failed")
    if max(abs(row["telescoping_balance_error"]) for row in energy_rows) >= 1e-10:
        raise AssertionError("a regularized polynomial-energy balance failed")


def main() -> None:
    check_open_cell_ode()
    check_event_scan()
    print("prime-power B-cell dynamics checks: PASS")


if __name__ == "__main__":
    main()
