#!/usr/bin/env python3
"""Fast regression tests for the finite E_h scaling-law audit."""

from __future__ import annotations

import csv
import json
import sys
import unittest
from pathlib import Path

import numpy as np

DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(DIRECTORY))

import audit_eh_scaling_law as audit


AUDIT_REPORT = DIRECTORY / "eh_scaling_law_audit_report.json"
ARCHIVED_SCAN = DIRECTORY / "exploratory_signed_triangular_scan_1e8_report.json"
SAMPLE_NAME = "signed_triangular_dynamics_lambda_0p25_sample.csv"
SAMPLE_CANDIDATES = (
    DIRECTORY.parent / "output" / SAMPLE_NAME,
    DIRECTORY.parent / "github_upload" / "papers" / "figures" / SAMPLE_NAME,
    DIRECTORY.parents[1] / "papers" / "figures" / SAMPLE_NAME,
)
PLOTTED_SAMPLE = next(
    (path for path in SAMPLE_CANDIDATES if path.is_file()),
    SAMPLE_CANDIDATES[0],
)


class ScalingAuditUnitTests(unittest.TestCase):
    """Check calculation helpers without repeating the full 10^8 scan."""

    def test_linear_fit_recovers_exact_slope(self) -> None:
        x = np.asarray([1.0, 2.0, 3.0, 4.0])
        y = 1.5 + 2.25 * x
        fit = audit.linear_fit(x, y)
        self.assertAlmostEqual(float(fit["slope"]), 2.25, places=14)
        self.assertAlmostEqual(float(fit["intercept"]), 1.5, places=14)
        self.assertAlmostEqual(float(fit["R_squared"]), 1.0, places=14)

    def test_quantile_record_has_expected_grain(self) -> None:
        record = audit.quantile_record(np.arange(1.0, 1001.0))
        self.assertEqual(record["count"], 1000)
        self.assertEqual(record["maximum"], 1000.0)
        self.assertAlmostEqual(record["p50"], 500.5)

    def test_saved_full_population_report_matches_archived_maxima(self) -> None:
        report = json.loads(AUDIT_REPORT.read_text(encoding="utf-8"))
        archived = json.loads(ARCHIVED_SCAN.read_text(encoding="utf-8"))
        archived_rows = archived["lambda_results"]["lambda_0p25"][
            "decade_diagnostics"
        ]
        archived_by_decade = {row["decade"]: row for row in archived_rows}

        self.assertEqual(report["population"]["center_count"], 5_762_859)
        self.assertEqual(report["population"]["prime_count"], 5_761_455)
        self.assertEqual(report["population"]["higher_prime_power_count"], 1_404)
        self.assertEqual(report["configuration"]["preferred_parsimonious_A"], 1.0)
        self.assertEqual(
            report["runtime"]["source_sha256"]["audit_eh_scaling_law.py"],
            audit.sha256(audit.SCRIPT_PATH),
        )

        for row in report["decades_for_A_equals_1"]:
            archived_maximum = archived_by_decade[row["decade"]]["extrema"][
                "maximum_abs_E_h_over_R"
            ]["value"]
            self.assertEqual(row["abs_E_over_Q"]["maximum"], archived_maximum)

    def test_saved_plot_sample_count_matches_csv(self) -> None:
        report = json.loads(AUDIT_REPORT.read_text(encoding="utf-8"))
        with PLOTTED_SAMPLE.open("r", encoding="utf-8", newline="") as handle:
            csv_count = sum(1 for _ in csv.DictReader(handle))
        self.assertEqual(
            report["configuration"]["plot_sample_actual_count"], csv_count
        )
        self.assertEqual(csv_count, 17_836)


if __name__ == "__main__":
    unittest.main()
