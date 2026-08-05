#!/usr/bin/env python3
"""Regression checks for the exact prime-side decomposition of E_h.

The tested formula is equation (3.3) in
Prime_Side_Eh_Decay_Derivation_2026-08-05.md.
"""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from exploratory_eventwise_reserve_scan import (  # noqa: E402
    EULER_GAMMA,
    LOG_TWO_PI,
    trivial_f_array,
)
from exploratory_signed_triangular_scan import (  # noqa: E402
    build_event_support,
    phi_at_x,
    required_support_limit,
)


def trivial_zero_remainder(q: int, h: float, terms: int = 80) -> float:
    """Return D_h R_0-R_0' in a numerically stable exponential form."""

    t = math.log(q)
    total = 0.0
    for k in range(1, terms + 1):
        alpha = 2.0 * k + 0.5
        coefficient = 1.0 / (2.0 * k * (2.0 * k + 1.0))
        centered_difference = (
            math.exp(-alpha * (t - h)) - math.exp(-alpha * (t + h))
        ) / (2.0 * h)
        total += coefficient * (
            alpha * math.exp(-alpha * t) - centered_difference
        )
    return total


class PrimeSideEhDecompositionTest(unittest.TestCase):
    """Compare equation (3.3) with independent endpoint evaluation."""

    centers = (2, 5, 16, 125, 997, 4096, 7919)
    lambdas = (0.25, 0.5, 1.0)

    @classmethod
    def setUpClass(cls) -> None:
        support_limit = required_support_limit(max(cls.centers), max(cls.lambdas))
        (
            cls.event_q,
            cls.event_prime,
            cls.event_exponent,
            cls.event_lambda,
        ) = build_event_support(support_limit, segment_span=20_000)
        cls.cumulative_psi = np.cumsum(cls.event_lambda, dtype=np.float64)
        cls.cumulative_s_lambda = np.cumsum(
            cls.event_lambda / cls.event_q.astype(np.float64),
            dtype=np.float64,
        )

    def evaluate_direct_and_decomposed(
        self, q: int, lambda_value: float
    ) -> tuple[float, float]:
        """Return E_h from endpoints and from the exact decomposition."""

        position = int(np.searchsorted(self.event_q, q, side="right") - 1)
        self.assertEqual(int(self.event_q[position]), q)

        t = math.log(q)
        root_q = math.sqrt(q)
        h = lambda_value * t / root_q
        lambda_q = float(self.event_lambda[position])
        psi_right = float(self.cumulative_psi[position])
        s_right = float(self.cumulative_s_lambda[position])

        reciprocal_q = 1.0 / q
        trivial_f = float(
            trivial_f_array(np.array([reciprocal_q], dtype=np.float64))[0]
        )
        trivial_f_prime = math.log1p(-(reciprocal_q * reciprocal_q))
        u_right = (
            s_right
            - t
            + EULER_GAMMA
            + 0.5 * reciprocal_q * trivial_f_prime
            - 0.5 * trivial_f
        )
        v_right = (
            psi_right - q + LOG_TWO_PI + 0.5 * trivial_f_prime
        )
        p_right = -0.5 * (root_q * u_right + v_right / root_q)
        p_bar = p_right + 0.5 * lambda_q / root_q

        endpoints = np.array(
            [q * math.exp(-h), q * math.exp(h)], dtype=np.float64
        )
        phi_minus, phi_plus = phi_at_x(
            endpoints,
            self.event_q,
            self.cumulative_psi,
            self.cumulative_s_lambda,
        )
        direct = p_bar - float(phi_plus - phi_minus) / (2.0 * h)

        s_star = s_right - 0.5 * lambda_q / q
        psi_star = psi_right - 0.5 * lambda_q
        u_star = s_star - t + EULER_GAMMA
        v_star = psi_star - q + LOG_TWO_PI
        b_hat_midpoint = 0.5 * (
            root_q * u_star + v_star / root_q
        )

        event_x = self.event_q.astype(np.float64)
        right_mask = (event_x > q) & (event_x < q * math.exp(h))
        left_mask = (event_x < q) & (event_x > q * math.exp(-h))

        right_x = event_x[right_mask]
        right_u = np.log(right_x / q)
        right_sum = np.sum(
            self.event_lambda[right_mask]
            * np.sqrt(q / right_x)
            * np.sinh((h - right_u) / 2.0),
            dtype=np.float64,
        )

        left_x = event_x[left_mask]
        left_u = np.log(q / left_x)
        left_sum = np.sum(
            self.event_lambda[left_mask]
            * np.sqrt(q / left_x)
            * np.sinh((h - left_u) / 2.0),
            dtype=np.float64,
        )
        delta_sinh = 2.0 * float(right_sum - left_sum) / h

        kappa = 2.0 * math.sinh(h / 2.0) / h - 1.0
        beta = (
            2.0 * math.sinh(h / 2.0) / h - math.cosh(h / 2.0)
        )
        decomposed = (
            delta_sinh / (2.0 * root_q)
            + kappa * b_hat_midpoint
            + root_q * beta
            + trivial_zero_remainder(q, h)
        )
        return direct, decomposed

    def test_exact_decomposition_at_prime_and_higher_power_centers(self) -> None:
        """Check several prime/higher-power centers and window multipliers."""

        for q in self.centers:
            for lambda_value in self.lambdas:
                with self.subTest(q=q, lambda_value=lambda_value):
                    direct, decomposed = self.evaluate_direct_and_decomposed(
                        q, lambda_value
                    )
                    # Prefix sums and endpoint subtraction are independent
                    # binary64 paths, so allow accumulated rounding at the
                    # largest center while remaining far below plot scale.
                    tolerance = 5.0e-11 * (1.0 + abs(direct))
                    self.assertLess(abs(direct - decomposed), tolerance)


if __name__ == "__main__":
    unittest.main()
