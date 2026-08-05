# Finite Robin verification and prime-power dynamics

This directory is the computational companion to
[`papers/Finite_Robin_Verification_via_CA_Prime_Power_Reduction_en.tex`](../../papers/Finite_Robin_Verification_via_CA_Prime_Power_Reduction_en.tex).

## Scope

The files have two separate roles:

- The finite-certificate programs use exact integer operations and directed
  `decimal` interval arithmetic. Their archived JSON reports support the two
  finite theorems in the paper.
- `exploratory_b_cell_dynamics.py` and the signed-triangular scanner investigate
  the signed prime-power dynamics in binary64 arithmetic. Their regression and
  independent direct-QA tests are exploratory and are not part of the finite
  interval certificate.
- `test_prime_side_eh_decomposition.py` checks the exact prime-side
  decomposition of the signed triangular discrepancy through independent
  endpoint and local-window evaluations. `audit_eh_scaling_law.py` then audits
  the proposed finite scaling on all 5,762,859 prime-power events through
  `10^8`. Its fitted exponents, quantiles, and archived JSON report are
  exploratory binary64 evidence, not an asymptotic theorem or interval
  certificate.
- `plot_signed_triangular_dynamics.py` reproduces the four-panel figure used in
  the paper from the archived signed-triangular report and a deterministic
  event sample.

## Requirements

- The standard-library finite certifiers support Python 3.10 or newer.
- The pinned exploratory and figure environment requires Python 3.11 or
  newer. The archived runs used CPython 3.13.5.
- The finite certificate programs use only the Python standard library.
- The exploratory programs additionally use NumPy. The archived scan used
  NumPy 2.3.5.
- The independent 80-digit signed-triangular QA additionally uses mpmath 1.3.0.
- Figure generation additionally uses Matplotlib 3.10.7.

Install the exploratory dependencies with:

```powershell
python -m pip install -r requirements.txt
```

## Non-destructive integrity audit

Run from the repository root:

```powershell
Set-Location verification/finite_robin_prime_power
python test_finite_robin_paper_manifest.py
python test_exploratory_b_cell_dynamics.py
python test_exploratory_signed_triangular_scan.py
python test_exploratory_signed_triangular_direct_qa.py
python test_prime_side_eh_decomposition.py
python test_audit_eh_scaling_law.py
```

The manifest validator reopens the archived reports, checks their status and
headline values, verifies cross-precision interval nesting, and compares the
SHA-256 hashes printed in the paper with the archived files. These commands do
not overwrite a canonical report. The final two tests respectively verify the
exact prime-side identity at prime and higher-prime-power centers and compare
the archived full-population scaling audit with the original signed scan and
the deterministic plotting CSV.

## Fresh recomputation without overwriting the archive

The commands below preserve the archived JSON files by writing new results to
`recomputed/`. The direct all-profile run, the full exploratory event scan, and
the full-population scaling audit are the most expensive commands.

```powershell
New-Item -ItemType Directory -Force recomputed | Out-Null

python test_robin_base_5041_55440_interval.py `
  --start 5041 --end 55440 --precision 44 --gamma-terms 200000 `
  --report recomputed/robin_base_5041_55440_interval_report.json

python test_ca_support_11_overlap_certificate.py `
  --precision 60 `
  --report recomputed/ca_support_11_overlap_certificate_report.json

python test_ca_all_profile_interval_certificate.py `
  --min-support 3299 --max-support 56048351 --precision 44 `
  --gamma-terms 200000 --series-cutoff 20000 `
  --higher-search-limit 100000 --summary-block-size 250000 `
  --progress-every 250000 --sieve-padding 1000 `
  --endpoint-crosscheck-report ca_endpoint_interval_certificate_report.json `
  --report recomputed/ca_all_profile_interval_certificate_report.json

python test_ca_exact_buffer_interval_certifier.py `
  --min-prime 3299 --max-prime 56048351 --target 0.7825 `
  --precision 40 --series-cutoff 1000 --fixed-scale-digits 70 `
  --direct-check-primes 5297,5303,10007 --progress-every 500000 `
  --report recomputed/ca_exact_buffer_interval_certificate_report.json

Copy-Item test_asymptotic_dstar_lower_bound.py recomputed/
python recomputed/test_asymptotic_dstar_lower_bound.py

python test_ca_residual_upward_extended_certificate.py `
  --precision 60 --gamma-terms 100000 `
  --report recomputed/ca_residual_upward_extended_certificate_report.json
python test_ca_residual_upward_extended_certificate.py `
  --precision 80 --gamma-terms 100000 `
  --report recomputed/ca_residual_upward_extended_certificate_p80_report.json

python test_ca_all_profile_interval_certificate.py `
  --min-support 11 --max-support 3299 --precision 44 `
  --gamma-terms 200000 --series-cutoff 20000 `
  --higher-search-limit 3000 --summary-block-size 250000 `
  --progress-every 0 --sieve-padding 1000 `
  --endpoint-crosscheck-report ca_endpoint_interval_certificate_report.json `
  --report recomputed/ca_all_profile_bridge_11_3299_report.json

python test_ca_all_integer_cutoff_certificate.py `
  --precision 60 `
  --report recomputed/ca_all_integer_cutoff_certificate_report.json
python test_ca_all_integer_cutoff_certificate.py `
  --precision 80 `
  --report recomputed/ca_all_integer_cutoff_certificate_p80_report.json

python exploratory_b_cell_dynamics.py `
  --limit 56048351 --scan-min 2 --bin-width 1 `
  --powers 0 1 2 3 4 6 8 `
  --report recomputed/exploratory_b_cell_dynamics_report.json

python exploratory_signed_triangular_scan.py `
  --limit 100000000 --segment-span 10000000 --chunk-size 250000 `
  --lambdas 0.25 0.30 0.50 1.00 `
  --report recomputed/exploratory_signed_triangular_scan_1e8_report.json

python audit_eh_scaling_law.py `
  --limit 100000000 --lambda-value 0.25 `
  --segment-span 10000000 --chunk-size 250000 --sample-count 30000 `
  --report recomputed/eh_scaling_law_audit_report.json
```

The dependent certifiers above intentionally validate the immutable archived
dependency chain. For a completely fresh end-to-end run, work in a disposable
copy of this directory and let each program write its canonical report name.

The JSON reports contain environment-dependent metadata, including absolute
paths, command-line arguments, platform strings, and timings. Consequently, a
fresh whole-file SHA-256 is not expected to equal the archived report hash even
when all certified mathematical fields agree.

## Reproduce the paper figure

From the repository root, run:

```powershell
Set-Location verification/finite_robin_prime_power
python plot_signed_triangular_dynamics.py `
  --limit 100000000 --lambda-value 0.25 --samples 30000 `
  --segment-span 10000000 `
  --report exploratory_signed_triangular_scan_1e8_report.json
```

The script writes the following files to `papers/figures/`:

- `signed_triangular_dynamics_lambda_0p25.png`, used by the TeX paper;
- `signed_triangular_dynamics_lambda_0p25.svg`, the vector version;
- `signed_triangular_dynamics_lambda_0p25_sample.csv`, the plotted sample.

The deterministic sample contains 17,836 prime-power events: logarithmically
spaced centers, every early event through 10,000, and every higher prime power
through the limit. The plotted decadal envelopes are not sampled; they are
read from the archived full scan over all 5,762,859 centers. The figure and
CSV are exploratory binary64 artifacts, not interval certificates.

## Included helper files

The flat directory also contains the transitive interval backends and the
endpoint cross-check report required by the public commands. The exploratory
prime-side helpers are included for the B-cell and signed-triangular
implementations, together with the paper-figure generator, the exact
decomposition test, the scaling audit and report, and their derivation and
interpretation notes. No zero table or external CSV input is required by the
commands above.
