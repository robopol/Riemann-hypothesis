# Robin-MVDC Status Support Files

This folder contains the auxiliary Python scripts and generated audit tables
used for the corrected Robin-MVDC status manuscript.

The files support the two open analytic routes described in that manuscript;
they are not presented as a proof.

The files are staged here for upload to:

https://github.com/robopol/Riemann-hypothesis

## Main script groups

- Signed Chebyshev / MVDC centre route:
  `test_w_mvdc_certified_envelope.py`,
  `test_w_mvdc_centered_m1_partition.py`,
  `test_w_chebyshev_target.py`,
  `test_mvdc_signed_chebyshev_center.py`,
  `test_mvdc_signed_chebyshev_ledger.py`,
  `test_w_mvdc_reduced_blocks.py`,
  `test_mvdc_prime_cell_moments.py`.

- CA-cap and beta-envelope route:
  `test_ca_cap_envelope.py`,
  `test_optimal_ledger_envelope.py`,
  `test_mvdc_beta_requirement.py`,
  `test_beta_env_gap_absorption.py`,
  `test_ca_like_ledger.py`.

- First-moment reciprocal-prime route:
  `test_ca_sminus_step_threshold.py`,
  `test_ca_cumulative_mvdc_envelope.py`,
  `test_ca_successive_mvdc_blocks.py`.

- Input/audit generators used by the scripts:
  `test_w_corrected_blocks.py`,
  `test_ca_theta_cancellation.py`,
  `test_mvdc_ca_beta_target.py`.

## Notes

The `.csv` and `.json` files are generated audit outputs kept beside the
scripts because several scripts use same-directory default input paths. Run
the scripts from this folder, or pass explicit `--input`, `--csv`, and `--json`
arguments if the files are moved.

These scripts are numerical audit tools. They check algebraic bookkeeping,
scales, and consistency of the proposed analytic routes; they are not a
replacement for the missing analytic theorem stated in the manuscript.
