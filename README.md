# Riemann Hypothesis: Robin-Inequality Research Archive

This repository preserves the earlier corrected Robin-MVDC status archive and
adds a new finite computer-assisted verification of Robin's inequality,
together with its source, reproducibility programs, interval certificates,
recorded reports, and exploratory prime-power diagnostics.

The new theorem is explicitly finite. Neither it nor the exploratory material
proves Robin's inequality for all integers or proves the Riemann Hypothesis.

## Latest finite verification (5 August 2026)

- [`papers/Finite_Robin_Verification_via_CA_Prime_Power_Reduction_en.pdf`](papers/Finite_Robin_Verification_via_CA_Prime_Power_Reduction_en.pdf)
  English paper: *A Finite Computer-Assisted Verification of Robin's
  Inequality via Colossally Abundant Profiles, with Exact Prime-Power Residual
  Dynamics*.

- [`papers/Finite_Robin_Verification_via_CA_Prime_Power_Reduction_en.tex`](papers/Finite_Robin_Verification_via_CA_Prime_Power_Reduction_en.tex)
  LaTeX source of the paper.

- [`papers/figures/`](papers/figures/)
  PNG and SVG versions of the signed-triangular dynamics figure, plus the
  deterministic plotted sample.

The paper gives an unconditional finite, computer-assisted verification of

```text
sigma(n) < exp(gamma) * n * log(log(n))
```

for every integer

```text
5041 <= n <= 10^(7.1e22).
```

The finite certificate combines exhaustive directed-interval enumeration of
3,341,978 colossally abundant exponent profiles, an analytic prime-power
reduction, a smoothed explicit formula, finite-height RH verification, and an
absolute bound for the remaining zero tail. The certified support calculation
extends through `1.64967e23`.

The later sections retain the signs discarded by the absolute estimate and
derive exact prime-side residual, cell, kick, and signed-triangular identities.
The associated scans through `10^8` are exploratory finite evidence. They do
not prove the still-missing uniform eventwise inequality.

## Earlier corrected-status papers

The following older files are retained as part of the research history:

- [`papers/RH_MVDC_corrected_status_2026_en.pdf`](papers/RH_MVDC_corrected_status_2026_en.pdf)
  English corrected-status manuscript. It withdraws the previous conclusive
  claim and records the remaining analytic targets at that stage.

- [`papers/RH_MVDC_corrected_status_2026_SK.pdf`](papers/RH_MVDC_corrected_status_2026_SK.pdf)
  Slovak version of the corrected-status manuscript.

- [`papers/Appendix_RH_en.pdf`](papers/Appendix_RH_en.pdf)
  Auxiliary logarithmic bridge appendix for `P(N) < log(N)` for a least
  hypothetical Robin counterexample `N > 5040`.

- [`papers/Derived_Prime_Harmonic_Envelope_on_CA_Support_en.pdf`](papers/Derived_Prime_Harmonic_Envelope_on_CA_Support_en.pdf)
  Earlier derived prime-harmonic-envelope paper retained for provenance.

Zenodo record for the earlier corrected-status package:

<https://zenodo.org/records/20329550>

## Verification material

- [`verification/finite_robin_prime_power/`](verification/finite_robin_prime_power/)
  Complete computational companion to the latest paper. It contains the finite
  interval certifiers, archived JSON reports, manifest/hash validator,
  exploratory prime-power programs, analytic investigation notes, regression
  tests, and figure generator. Its own README separates rigorous finite
  certificates from exploratory binary64 diagnostics.

- [`verification/verify_low_prime_robin_region.py`](verification/verify_low_prime_robin_region.py)
  Earlier finite low-prime-region audit.

- [`verification/robin_mvdc_status_support/`](verification/robin_mvdc_status_support/)
  Scripts and recorded audit data for the corrected Robin-MVDC status routes.

- [`verification/strongest_theorem_support/`](verification/strongest_theorem_support/)
  Reproducibility files for strongest-beta and structured-candidate numerical
  checks.

The recorded JSON, CSV, and figure files document the archived runs. Numerical
or exploratory scripts do not replace the missing analytic theorem needed for
an infinite conclusion.

## Reproduction

Run the non-destructive integrity audit for the latest paper from the
repository root:

```powershell
Set-Location verification/finite_robin_prime_power
python test_finite_robin_paper_manifest.py
```

Install the optional exploratory and plotting dependencies with:

```powershell
python -m pip install -r requirements.txt
```

Additional tests and fresh recomputation commands are documented in
[`verification/finite_robin_prime_power/README.md`](verification/finite_robin_prime_power/README.md).
The expensive full-profile and full-event scans are not required for the
non-destructive manifest audit.

Earlier low-prime audit:

```powershell
Set-Location ../..
python verification/verify_low_prime_robin_region.py
```

Selected corrected-status audits:

```powershell
Set-Location verification/robin_mvdc_status_support
python test_mvdc_signed_chebyshev_center.py
python test_ca_sminus_step_threshold.py
```

Selected strongest-theorem support checks:

```powershell
Set-Location ../strongest_theorem_support
python scripts/run_sigma_max_beta_test.py --max-count 1000000 --checkpoints 24
python scripts/run_hcn_sa_bridge_test.py --max-count 40000 --checkpoints 10
```
