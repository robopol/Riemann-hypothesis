# Riemann Hypothesis: Robin-MVDC Corrected Status Archive

This repository contains the cleaned public support archive for the corrected
Robin-MVDC status manuscript. It is presented as a status archive with open
analytic targets, not as a closed solution of Robin's inequality or the Riemann
Hypothesis.

The earlier Robin-Mertens compensation branch and the consolidated proof text
are superseded. The current package keeps the valid logarithmic bridge appendix,
the corrected status manuscript, and numerical audit scripts/data for the two
remaining analytic routes.

## Current Papers

- `papers/RH_MVDC_corrected_status_2026_en.pdf`  
  Corrected status manuscript. It withdraws the previous conclusive claim and
  reformulates the remaining work as two analytic targets: a signed
  Chebyshev/MVDC centre route and a first-moment reciprocal-prime certificate.

- `papers/RH_MVDC_corrected_status_2026_SK.pdf`  
  Slovak version of the corrected status manuscript.

- `papers/Appendix_RH_en.pdf`  
  Auxiliary logarithmic bridge appendix for
  `P(N) < log N` for a least hypothetical Robin counterexample `N > 5040`.

Zenodo record:

<https://zenodo.org/records/20329550>

## Verification Material

- `verification/verify_low_prime_robin_region.py`  
  Finite low-prime region audit used to check the explicit small-prime range.

- `verification/robin_mvdc_status_support/`  
  Scripts and recorded audit data for the corrected Robin-MVDC status routes.

- `verification/strongest_theorem_support/`  
  Reproducibility files for strongest-beta and structured-candidate numerical
  checks.

These scripts check algebraic bookkeeping, numerical scale, and
reproducibility. They do not replace the missing analytic theorem stated in the
corrected status manuscript.

## Reproduction

Run from the repository root:

```powershell
py verification/verify_low_prime_robin_region.py
```

Selected Robin-MVDC status audits:

```powershell
cd verification/robin_mvdc_status_support
py test_mvdc_signed_chebyshev_center.py
py test_ca_sminus_step_threshold.py
```

Selected strongest-theorem support checks:

```powershell
cd ../strongest_theorem_support
py scripts/run_sigma_max_beta_test.py --max-count 1000000 --checkpoints 24
py scripts/run_hcn_sa_bridge_test.py --max-count 40000 --checkpoints 10
```

The JSON, CSV, and figure files in the verification folders are recorded outputs
from these numerical runs.
