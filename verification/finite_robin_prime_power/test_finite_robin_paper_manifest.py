#!/usr/bin/env python3
"""Validate the standalone finite Robin paper against its certificate artifacts."""

from __future__ import annotations

import hashlib
import json
import re
import sys
from decimal import Decimal
from pathlib import Path
from typing import Any


SCRIPT_PATH = Path(__file__).resolve()
REPOSITORY_ROOT = SCRIPT_PATH.parents[2]
SUPPORT_ROOT = SCRIPT_PATH.parent
PUBLIC_ARTIFACT_ROOT = "verification/finite_robin_prime_power"
PAPER_NAME = "Finite_Robin_Verification_via_CA_Prime_Power_Reduction_en.tex"
PAPER_CANDIDATES = (
    REPOSITORY_ROOT / "papers" / PAPER_NAME,
    REPOSITORY_ROOT / "current_papers" / PAPER_NAME,
)
PAPER_PATH = next((path for path in PAPER_CANDIDATES if path.is_file()), PAPER_CANDIDATES[0])

EXPECTED_HASHES = {
    "current_support/robin_mvdc_status_support/test_ca_all_profile_interval_certificate.py":
        "8da8ad3d63b6dec8ffd4838b968a0added4c26eb5810fa2969dadb917f691077",
    "current_support/robin_mvdc_status_support/ca_all_profile_interval_certificate_report.json":
        "8ac1ee276df64cc91d93b483b501d34f8a9b69408d2fb0650d4bab8f43304f83",
    "current_support/robin_mvdc_status_support/test_ca_exact_buffer_interval_certifier.py":
        "add8f04a2e02b606491c1bdf6aab4fcb23f5d10e2db202326fded0c5035d5ef2",
    "current_support/robin_mvdc_status_support/ca_exact_buffer_interval_certificate_report.json":
        "57c14d33aa47b7ddff012043304151fc9b11b50a114c2634665f7455ce80dba1",
    "current_support/robin_mvdc_status_support/test_asymptotic_dstar_lower_bound.py":
        "c7981a2fdec34211e07434d82931c2ac534a20e384fdace31b4ee6d79a8211b5",
    "current_support/robin_mvdc_status_support/asymptotic_dstar_lower_bound_report.json":
        "ea33332d7acade20b35409131a6a5acac00d722cf3439d9b36cb7992d296fa18",
    "current_support/robin_mvdc_status_support/test_ca_residual_upward_extended_certificate.py":
        "165f7a069c93ede79709ed559fdd1a9908e41230fd4b520844931789a97d3e95",
    "current_support/robin_mvdc_status_support/ca_residual_upward_extended_certificate_report.json":
        "938713a1e1443196e349269d9f88f429d3859aa262fbf11042733c6ef1278263",
    "current_support/robin_mvdc_status_support/ca_residual_upward_extended_certificate_p80_report.json":
        "41244347c6366c3321dc653235b9db32008e38678654d2f899838b7f5c2ee8a1",
    "current_support/robin_mvdc_status_support/test_robin_base_5041_55440_interval.py":
        "ee296aa4a9e57e576da69a114e5983d9dad172d89e35d1498c83a1cff211c1aa",
    "current_support/robin_mvdc_status_support/robin_base_5041_55440_interval_report.json":
        "161381cb69f55030f2f6bb0c57c6518383aa66a7d65165cf98e1686aa013cf45",
    "current_support/robin_mvdc_status_support/test_ca_support_11_overlap_certificate.py":
        "660420eba61fa0ad61df16ed6117dfd76ac84adb78cff2892c44e0b500de29a9",
    "current_support/robin_mvdc_status_support/ca_support_11_overlap_certificate_report.json":
        "af903d574499c2b9450b9afed9b4799d3d01a668eb30097306d61ba97da1a0bd",
    "current_support/robin_mvdc_status_support/ca_all_profile_bridge_11_3299_report.json":
        "e29d66dde12af0e773d86827b0237ed323f62abfd19bc7fe028a7e065782e5ca",
    "current_support/robin_mvdc_status_support/test_ca_all_integer_cutoff_certificate.py":
        "df1238bee0184a3dccebdcaf72920b9362984ec845273692a786bf69030e4052",
    "current_support/robin_mvdc_status_support/ca_all_integer_cutoff_certificate_report.json":
        "0ba63f2645fc5c312aa00890daff8b87b720e68656c5cddfc31dd1465284f8f3",
    "current_support/robin_mvdc_status_support/ca_all_integer_cutoff_certificate_p80_report.json":
        "5dde972048d365efe11a9305bec698a0f0c56de0e27957f8d72f67a72b00351e",
    "current_support/robin_mvdc_status_support/exploratory_b_cell_dynamics.py":
        "1f1722248765f2ba314fbb5c17cd7727a58ef363f00254ed9ff9edb0fbd68c21",
    "current_support/robin_mvdc_status_support/exploratory_b_cell_dynamics_report.json":
        "5c610b2e08ccc94b30b158d1cd578dafad3189042491aa4e6d7006670c25da41",
    "current_support/robin_mvdc_status_support/test_exploratory_b_cell_dynamics.py":
        "446f0e22fe617c1ed7809dd3d2e78f1a10420cd611be4c8b33239ceb5e17b773",
    "working_notes/B_Prime_Power_Cell_Dynamics_Investigation_2026-08-04.md":
        "acb81a8c38101b3e5fada926faf75ac24a05e67d7d2315667c6aac49bff22b8b",
    "verification/finite_robin_prime_power/audit_eh_scaling_law.py":
        "f8c3d03963705cddead22ab07ff67a8f8d17855ee2d363cb8d0f5d0eebe6100b",
    "verification/finite_robin_prime_power/test_audit_eh_scaling_law.py":
        "f2987d0407dd59af4d0606d12f617ba50ff7e0afb9f587d78c12116226e35d6d",
    "verification/finite_robin_prime_power/eh_scaling_law_audit_report.json":
        "0aaf4eb5dd456e94417853dcf4609d43c0113fb5b14a097928589e9f56d3fa4b",
    "verification/finite_robin_prime_power/test_prime_side_eh_decomposition.py":
        "6ba6b5e18067cc2416df8624341a29d4c5b4d47910ce1f5e21c264186d4837bd",
    "verification/finite_robin_prime_power/Eh_Decay_Analytic_Investigation_2026-08-05.md":
        "de47320a9142d0142f6dbdfdab60e6eae6c1a935f923530a69c39d4166fee4c1",
    "verification/finite_robin_prime_power/Prime_Side_Eh_Decay_Derivation_2026-08-05.md":
        "2a05e34cef9ae3b0b5693e70b0da6eb2f33d118cd56385bdb6708616fe1b38b8",
    "verification/finite_robin_prime_power/Eh_Scaling_Law_Audit_2026-08-05.md":
        "11132a22a95e0937fc323739b3ef63e94febb9455e14f1e94d616b8f9fec77a8",
    "verification/finite_robin_prime_power/Spectral_Eh_Decay_Analysis_2026-08-05.md":
        "35c2f9a73b21969af7c9e893e7dc2f54170926b74aaff63bba0de165770f5870",
}

EXPECTED_MACROS = {
    "Xmin": "3299",
    "Xdirect": "56048351",
    "Xmax": "164967000000000000000000",
    "Kcut": "71000000000000000000000",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_report(name: str) -> dict[str, Any]:
    with (SUPPORT_ROOT / name).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def collect_intervals(value: Any, path: tuple[str, ...] = ()) -> dict[tuple[str, ...], tuple[Decimal, Decimal]]:
    intervals: dict[tuple[str, ...], tuple[Decimal, Decimal]] = {}
    if isinstance(value, dict):
        if "lower" in value and "upper" in value:
            try:
                intervals[path] = (Decimal(str(value["lower"])), Decimal(str(value["upper"])))
            except Exception:
                pass
        for key, child in value.items():
            intervals.update(collect_intervals(child, path + (str(key),)))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            intervals.update(collect_intervals(child, path + (str(index),)))
    return intervals


def main() -> int:
    failures: list[str] = []
    check_count = 0

    def require(condition: bool, message: str) -> None:
        nonlocal check_count
        check_count += 1
        if not condition:
            failures.append(message)

    paper = PAPER_PATH.read_text(encoding="utf-8")
    paper_without_verbatim = re.sub(
        r"\\begin\{verbatim\}.*?\\end\{verbatim\}",
        "",
        paper,
        flags=re.DOTALL,
    )

    environment_stack: list[str] = []
    environment_error: str | None = None
    for action, environment in re.findall(r"\\(begin|end)\{([^}]+)\}", paper_without_verbatim):
        if action == "begin":
            environment_stack.append(environment)
        elif not environment_stack or environment_stack[-1] != environment:
            environment_error = (
                f"Environment nesting mismatch at end of {environment}; "
                f"stack is {environment_stack}"
            )
            break
        else:
            environment_stack.pop()
    require(environment_error is None, environment_error or "Environment nesting is invalid")
    require(not environment_stack, f"Unclosed TeX environments: {environment_stack}")

    brace_depth = 0
    minimum_brace_depth = 0
    for index, character in enumerate(paper_without_verbatim):
        if character not in "{}":
            continue
        preceding_backslashes = 0
        cursor = index - 1
        while cursor >= 0 and paper_without_verbatim[cursor] == "\\":
            preceding_backslashes += 1
            cursor -= 1
        if preceding_backslashes % 2:
            continue
        brace_depth += 1 if character == "{" else -1
        minimum_brace_depth = min(minimum_brace_depth, brace_depth)
    require(
        brace_depth == 0 and minimum_brace_depth == 0,
        f"TeX grouping-brace imbalance: depth={brace_depth}, minimum={minimum_brace_depth}",
    )

    labels = re.findall(r"\\label\{([^}]+)\}", paper_without_verbatim)
    references = set(re.findall(r"\\(?:eqref|ref)\{([^}]+)\}", paper_without_verbatim))
    require(len(labels) == len(set(labels)), "The paper contains duplicate labels")
    require(not (references - set(labels)), f"Missing TeX labels: {sorted(references - set(labels))}")

    bibliography_keys = set(re.findall(r"\\bibitem\{([^}]+)\}", paper_without_verbatim))
    citation_keys: set[str] = set()
    for citation_group in re.findall(r"\\cite(?:\[[^]]*\])?\{([^}]+)\}", paper_without_verbatim):
        citation_keys.update(key.strip() for key in citation_group.split(","))
    require(not (citation_keys - bibliography_keys), f"Missing bibliography entries: {sorted(citation_keys - bibliography_keys)}")
    require(not (bibliography_keys - citation_keys), f"Uncited bibliography entries: {sorted(bibliography_keys - citation_keys)}")
    require(chr(96) not in paper, "The TeX source contains a literal backtick")

    paper_hashes = {
        "".join(groups)
        for groups in re.findall(
            r"\\sha\{([0-9a-f]{16})\}\{([0-9a-f]{16})\}\{([0-9a-f]{16})\}\{([0-9a-f]{16})\}",
            paper,
        )
    }
    paper_paths = set(re.findall(r"\\path\{([^}]+)\}", paper))

    for macro, expected_digits in EXPECTED_MACROS.items():
        match = re.search(rf"\\newcommand\{{\\{macro}\}}\{{([^}}]+)\}}", paper)
        require(match is not None, f"Missing TeX macro {macro}")
        if match is not None:
            observed_digits = re.sub(r"\D", "", match.group(1))
            require(observed_digits == expected_digits, f"TeX macro {macro} has {observed_digits}, expected {expected_digits}")

    for relative_path, expected_hash in EXPECTED_HASHES.items():
        artifact_name = Path(relative_path).name
        artifact_path = SUPPORT_ROOT / artifact_name
        if not artifact_path.is_file():
            artifact_path = REPOSITORY_ROOT / Path(relative_path)
        public_path = f"{PUBLIC_ARTIFACT_ROOT}/{artifact_name}"
        require(artifact_path.is_file(), f"Missing artifact: {public_path}")
        if artifact_path.is_file():
            observed_hash = sha256(artifact_path)
            require(observed_hash == expected_hash, f"SHA-256 mismatch for {public_path}: {observed_hash}")
        require(public_path in paper_paths, f"Paper does not list artifact path: {public_path}")
        require(expected_hash in paper_hashes, f"Paper does not list artifact hash: {expected_hash}")

    own_public_path = f"{PUBLIC_ARTIFACT_ROOT}/{SCRIPT_PATH.name}"
    own_hash = sha256(SCRIPT_PATH)
    require(own_public_path in paper_paths, "Paper does not list this manifest validator")
    require(own_hash in paper_hashes, f"Paper does not list current manifest-validator hash: {own_hash}")

    prime_stream_hash = "9438df02de33467c1ed4307c0a801c2c57777e00bb7413b17ce955a117c8b2c9"
    require(prime_stream_hash in paper_hashes, "Paper does not list the common prime-stream hash")

    direct = load_report("ca_all_profile_interval_certificate_report.json")
    require(direct["certificate_status"] == "PASS", "Direct CA report is not PASS")
    require(direct["configuration"]["min_support_prime"] == 3299, "Direct CA lower support mismatch")
    require(direct["configuration"]["max_support_prime"] == 56048351, "Direct CA upper support mismatch")
    require(direct["profile_coverage"]["prime_support_count"] == 3340551, "Direct support count mismatch")
    require(direct["profile_coverage"]["evaluated_distinct_numerical_profile_count"] == 3341978, "Direct profile count mismatch")
    require(direct["classification"]["active_higher_transition_ties"] == 1427, "Direct tie count mismatch")
    require(direct["profile_coverage"]["nonpositive_or_unresolved_profile_count"] == 0, "Direct CA report has a failed or unresolved profile")
    require(direct["prime_stream_sha256_little_endian_uint64_through_max_support"] == prime_stream_hash, "Direct prime-stream hash mismatch")

    dstar = load_report("ca_exact_buffer_interval_certificate_report.json")
    require(dstar["status"] == "PASS", "Finite D-star report is not PASS")
    require(dstar["status_counts"] == {"support_count": 3340551, "pass": 3340551, "fail": 0, "inconclusive": 0}, "Finite D-star status counts mismatch")
    require(Decimal(dstar["minimum"]["scaled_buffer"]["lower"]) > Decimal("0.7825"), "Finite D-star minimum is too small")

    asymptotic = load_report("asymptotic_dstar_lower_bound_report.json")
    require(asymptotic["status"] == "PASS", "Asymptotic D-star report is not PASS")
    require(Decimal(asymptotic["endpoint"]["scaled_dstar_lower"]["lower"]) > Decimal("0.5161"), "Asymptotic D-star endpoint is too small")

    residual_60 = load_report("ca_residual_upward_extended_certificate_report.json")
    residual_80 = load_report("ca_residual_upward_extended_certificate_p80_report.json")
    require(residual_60["status"] == residual_80["status"] == "PASS", "Extended residual report is not PASS at both precisions")
    require(residual_60["coverage"]["segment_count"] == 662, "Extended residual segment count mismatch")
    minimum_segment = residual_60["coverage"]["minimum_segment"]
    require(minimum_segment["left"] == 164950000000000000000000, "Minimum segment left endpoint mismatch")
    require(minimum_segment["right"] == 164967000000000000000000, "Minimum segment right endpoint mismatch")
    require(Decimal(minimum_segment["combined_clearance"]["lower"]) > Decimal("2.1540734262e-7"), "Extended residual clearance is too small")
    barrier = residual_60["scalar_envelope_barrier_audit"]
    require(barrier["positive_point"] == 164967000000000000000000, "Positive barrier point mismatch")
    require(barrier["negative_point"] == 164968000000000000000000, "Negative barrier point mismatch")

    base = load_report("robin_base_5041_55440_interval_report.json")
    require(base["certificate_status"] == "PASS", "Integer base report is not PASS")
    require(base["coverage"]["integer_count"] == 50400 and base["coverage"]["failure_count"] == 0, "Integer base coverage mismatch")
    require(base["minimum_gap"]["n"] == 10080, "Integer base minimum location mismatch")

    overlap = load_report("ca_support_11_overlap_certificate_report.json")
    require(overlap["certificate_status"] == "PASS", "Support-11 overlap report is not PASS")
    require(overlap["old_tie_profile"]["integer"] == 5040, "Old support-11 tie profile mismatch")
    require(overlap["forced_exact_support_11_profile"]["integer"] == 55440, "Forced support-11 profile mismatch")
    require(not overlap["unresolved_comparisons"], "Support-11 overlap has unresolved comparisons")

    bridge = load_report("ca_all_profile_bridge_11_3299_report.json")
    require(bridge["certificate_status"] == "PASS", "Small CA bridge report is not PASS")
    require(bridge["profile_coverage"]["prime_support_count"] == 459, "Small CA bridge support count mismatch")
    require(bridge["profile_coverage"]["evaluated_distinct_numerical_profile_count"] == 503, "Small CA bridge profile count mismatch")
    require(bridge["classification"]["active_higher_transition_ties"] == 44, "Small CA bridge tie count mismatch")

    cutoff_60 = load_report("ca_all_integer_cutoff_certificate_report.json")
    cutoff_80 = load_report("ca_all_integer_cutoff_certificate_p80_report.json")
    require(cutoff_60["certificate_status"] == cutoff_80["certificate_status"] == "PASS", "All-integer cutoff report is not PASS at both precisions")
    require(cutoff_60["numeric_cutoff"]["X"] == EXPECTED_MACROS["Xmax"], "All-integer cutoff X mismatch")
    require(cutoff_60["numeric_cutoff"]["clean_decimal_exponent_K"] == EXPECTED_MACROS["Kcut"], "All-integer exponent K mismatch")

    exploratory = load_report("exploratory_b_cell_dynamics_report.json")
    require(exploratory["checks"]["all_pass"], "Exploratory B-cell report has a failed consistency check")
    require(exploratory["input"]["scan_limit"] == 56048351, "Exploratory B-cell scan limit mismatch")
    require(exploratory["input"]["all_prime_power_event_count"] == 3342115, "Exploratory B-cell event count mismatch")
    require(exploratory["stationary_extrema"]["stationary_cell_count"] == 23746, "Exploratory stationary-cell count mismatch")
    cone = exploratory["Phi_first_crossing_energy"]["negative_H_cone_margins"]
    require(cone["both_one_sided_states_strictly_inside_cone"], "Exploratory cone diagnostic failed")
    require(cone["minimum_pre_kick_ratio_event"] == 5, "Exploratory minimum kick-ratio event mismatch")
    require(
        Decimal(str(cone["minimum_pre_kick_ratio_R_plus_Phi_prime_left_over_j"]))
        > Decimal("1.2055"),
        "Exploratory minimum kick ratio is too small",
    )
    require(cone["minimum_event_lower_boundary_margin_event"] == 2, "Exploratory minimum post-kick reserve event mismatch")
    require(
        Decimal(str(cone["minimum_event_lower_boundary_margin_R_plus_Phi_prime_right"]))
        > Decimal("0.1212"),
        "Exploratory minimum post-kick reserve is too small",
    )
    regularized_rows = exploratory["regularized_Bhat_polynomial_energy"]["rows"]
    require([row["r"] for row in regularized_rows] == [0, 1, 2, 3, 4, 6, 8], "Exploratory polynomial-energy powers mismatch")
    require(
        all(row["signed_jump_plus_cell_sum"] < 0 for row in regularized_rows if row["r"] > 0),
        "An exploratory positive-power energy has nonnegative total signed change",
    )
    require(
        exploratory["runtime"]["script_sha256"]
        == EXPECTED_HASHES["current_support/robin_mvdc_status_support/exploratory_b_cell_dynamics.py"],
        "Exploratory report runtime script hash mismatch",
    )

    scaling = load_report("eh_scaling_law_audit_report.json")
    scaling_configuration = scaling["configuration"]
    scaling_population = scaling["population"]
    require(
        scaling["status"].startswith("Exploratory binary64 finite-range scaling audit."),
        "Scaling report status does not identify an exploratory binary64 audit",
    )
    require(scaling_configuration["limit"] == 100000000, "Scaling audit limit mismatch")
    require(Decimal(str(scaling_configuration["lambda"])) == Decimal("0.25"), "Scaling audit lambda mismatch")
    require(scaling_configuration["plot_sample_actual_count"] == 17836, "Scaling plot-sample count mismatch")
    require(Decimal(str(scaling_configuration["preferred_parsimonious_A"])) == Decimal("1"), "Scaling preferred logarithmic exponent mismatch")
    require(
        Decimal("1.12")
        < Decimal(str(scaling_configuration["median_fit_A_q_ge_1e4"]))
        < Decimal("1.14"),
        "Scaling median fit at q>=1e4 is outside its archived range",
    )
    require(
        Decimal("0.96")
        < Decimal(str(scaling_configuration["median_fit_A_q_ge_1e5"]))
        < Decimal("0.98"),
        "Scaling median fit at q>=1e5 is outside its archived range",
    )

    expected_scaling_population = {
        "center_count": 5762859,
        "prime_count": 5761455,
        "higher_prime_power_count": 1404,
        "positive_E_count": 2890468,
        "negative_E_count": 2872391,
        "zero_E_binary64_count": 0,
    }
    for field, expected_value in expected_scaling_population.items():
        require(
            scaling_population[field] == expected_value,
            f"Scaling population field {field} is {scaling_population[field]}, expected {expected_value}",
        )
    require(
        scaling_population["prime_count"] + scaling_population["higher_prime_power_count"]
        == scaling_population["center_count"],
        "Scaling prime and higher-power counts do not partition the event population",
    )
    require(
        scaling_population["positive_E_count"]
        + scaling_population["negative_E_count"]
        + scaling_population["zero_E_binary64_count"]
        == scaling_population["center_count"],
        "Scaling sign counts do not partition the event population",
    )

    scaling_script_key = "verification/finite_robin_prime_power/audit_eh_scaling_law.py"
    scaling_source_hashes = scaling["runtime"]["source_sha256"]
    require(
        scaling_source_hashes["audit_eh_scaling_law.py"] == EXPECTED_HASHES[scaling_script_key],
        "Scaling report audit-script hash mismatch",
    )
    for dependency_name in (
        "exploratory_signed_triangular_scan.py",
        "exploratory_eventwise_reserve_scan.py",
        "plot_signed_triangular_dynamics.py",
    ):
        dependency_path = SUPPORT_ROOT / dependency_name
        require(dependency_path.is_file(), f"Missing scaling dependency: {dependency_name}")
        if dependency_path.is_file():
            require(
                scaling_source_hashes[dependency_name] == sha256(dependency_path),
                f"Scaling report dependency hash mismatch for {dependency_name}",
            )

    signed_scan = load_report("exploratory_signed_triangular_scan_1e8_report.json")
    archived_decades = signed_scan["lambda_results"]["lambda_0p25"]["decade_diagnostics"]
    scaling_decades = scaling["decades_for_A_equals_1"]
    archived_by_decade = {row["decade"]: row for row in archived_decades}
    scaling_by_decade = {row["decade"]: row for row in scaling_decades}
    require(len(archived_decades) == len(scaling_decades) == 8, "Scaling decade-row count mismatch")
    require(
        set(archived_by_decade) == set(scaling_by_decade),
        "Scaling and signed-scan reports cover different decades",
    )
    for decade in sorted(set(archived_by_decade) & set(scaling_by_decade)):
        archived_maximum = archived_by_decade[decade]["extrema"]["maximum_abs_E_h_over_R"]["value"]
        scaling_row = scaling_by_decade[decade]
        scaling_maximum = scaling_row["abs_E_over_Q"]["maximum"]
        require(
            scaling_maximum == archived_maximum,
            f"Scaling maximum differs from the signed scan in decade {decade}",
        )
        normalized = scaling_row["Z_A_full"]
        require(
            normalized["count"] == scaling_row["event_count"],
            f"Scaling normalized count mismatch in decade {decade}",
        )
        require(
            0 < normalized["p50"] <= normalized["p90"] <= normalized["p99"]
            <= normalized["p99.9"] <= normalized["maximum"],
            f"Scaling normalized quantiles are not ordered in decade {decade}",
        )
    require(
        sum(row["event_count"] for row in scaling_decades) == scaling_population["center_count"],
        "Scaling decade counts do not sum to the full event population",
    )

    for label, low_report, high_report, expected_count in (
        ("extended residual", residual_60, residual_80, 4007),
        ("all-integer cutoff", cutoff_60, cutoff_80, 7),
    ):
        low_intervals = collect_intervals(low_report)
        high_intervals = collect_intervals(high_report)
        require(set(low_intervals) == set(high_intervals), f"{label} precision reports have different interval paths")
        require(len(low_intervals) == expected_count, f"{label} interval count is {len(low_intervals)}, expected {expected_count}")
        for interval_path in set(low_intervals) & set(high_intervals):
            low_lower, low_upper = low_intervals[interval_path]
            high_lower, high_upper = high_intervals[interval_path]
            require(
                low_lower <= high_lower <= high_upper <= low_upper,
                f"{label} high-precision interval is not nested at {'/'.join(interval_path)}",
            )

    if failures:
        print(f"FAIL: {len(failures)} of {check_count} checks failed")
        for failure in failures:
            print(f"- {failure}")
        return 1

    print(f"PASS: {check_count} paper-manifest, report, hash, and nesting checks")
    print("Extended residual intervals nested: 4007/4007")
    print("All-integer cutoff intervals nested: 7/7")
    print(f"Manifest validator SHA-256: {own_hash}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
