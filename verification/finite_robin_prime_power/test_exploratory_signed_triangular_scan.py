#!/usr/bin/env python3
"""Regression tests for the signed triangular event scan."""

from __future__ import annotations

import argparse
import math
import sys
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from exploratory_signed_triangular_scan import (  # noqa: E402
    EULER_GAMMA,
    LOG_TWO_PI,
    WEIL_C,
    required_support_limit,
    scan,
)


class SignedTriangularScanTests(unittest.TestCase):
    """Check finite counts, identities, and early limiting events."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.report = scan(
            argparse.Namespace(
                limit=100_000,
                segment_span=20_000,
                chunk_size=20_000,
                lambdas=(0.25, 0.30, 0.50, 1.00),
                progress_every=0,
            )
        )

    def test_event_count_and_internal_checks(self) -> None:
        self.assertEqual(
            self.report["counts"]["center_all_prime_power_events"], 9_700
        )
        self.assertTrue(self.report["checks"]["all_pass"])

    def test_turan_defect_is_positive_on_small_scan(self) -> None:
        for result in self.report["lambda_results"].values():
            self.assertEqual(result["raw_negative_U_count"], 0)
            self.assertGreater(
                result["global_extrema"]["minimum_U"]["value"], 0.0
            )

    def test_early_nonpositive_conditional_reserve_counts(self) -> None:
        expected = {
            "lambda_0p25": (7, 19),
            "lambda_0p3": (5, 19),
            "lambda_0p5": (3, 5),
            "lambda_1": (0, None),
        }
        for key, (negative_count, last_q) in expected.items():
            result = self.report["lambda_results"][key]
            self.assertEqual(
                result["M_T_sign_counts"]["negative"], negative_count
            )
            witness = result["global_extrema"].get(
                "maximum_q_with_nonpositive_M_T"
            )
            if last_q is None:
                self.assertIsNone(witness)
            else:
                self.assertEqual(witness["row"]["q"], last_q)

    def test_exact_endpoint_identities_close(self) -> None:
        for result in self.report["lambda_results"].values():
            self.assertLess(
                result["maximum_algebra_identity_residual"], 1e-12
            )
            self.assertLess(
                result["maximum_post_slope_identity_residual"], 1e-14
            )

    def test_lambda_one_is_positive_from_q_two(self) -> None:
        result = self.report["lambda_results"]["lambda_1"]
        minimum = result["global_extrema"]["minimum_M_T"]
        self.assertEqual(minimum["row"]["q"], 2)
        self.assertGreater(minimum["value"], 0.025)

    def test_phi_2h_against_explicit_no_prime_formula(self) -> None:
        """Check Phi near zero without using the prefix lookup evaluator."""

        row = self.report["lambda_results"]["lambda_0p25"][
            "global_extrema"
        ]["minimum_U"]["row"]
        s = 2.0 * row["h"]
        x = math.exp(s)
        self.assertLess(x, 2.0)
        root_x = math.sqrt(x)
        y = 1.0 / x
        f_value = (
            (1.0 + y) * math.log1p(y)
            + (y - 1.0) * math.log1p(-y)
            - 2.0 * y
        )
        b_no_prime = (
            root_x * (EULER_GAMMA - s + 1.0)
            - LOG_TWO_PI / root_x
            - 0.5 * root_x * f_value
        )
        phi_no_prime = WEIL_C - b_no_prime
        self.assertAlmostEqual(phi_no_prime, row["Phi_2h"], places=14)

    def test_support_bound_covers_large_window_multiplier(self) -> None:
        """Cover both endpoint families outside the production lambda grid."""

        center_limit = 100
        lambda_value = 10.0
        support_limit = required_support_limit(center_limit, lambda_value)
        for q in range(2, center_limit + 1):
            t = math.log(q)
            h = lambda_value * t / math.sqrt(q)
            self.assertLessEqual(q * math.exp(h), support_limit)
            self.assertLessEqual(math.exp(2.0 * h), support_limit)

    def test_scan_rejects_windows_that_cross_t_zero(self) -> None:
        """Require t-h>0 at the smallest center q=2."""

        with self.assertRaisesRegex(ValueError, "0<lambda<sqrt\\(2\\)"):
            scan(
                argparse.Namespace(
                    limit=100,
                    segment_span=100,
                    chunk_size=100,
                    lambdas=(math.sqrt(2.0),),
                    progress_every=0,
                )
            )


if __name__ == "__main__":
    unittest.main()
