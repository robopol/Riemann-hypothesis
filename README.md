# Riemann Hypothesis: Robin-Inequality Research Archive

This repository preserves the earlier corrected Robin-MVDC status archive and
adds a new finite computer-assisted verification of Robin's inequality,
together with its source, reproducibility programs, interval certificates,
recorded reports, and exploratory prime-power diagnostics.

The August study proves Robin's inequality on a finite range. The two November
preprints develop a prime-power criterion for RH and connect it to bounds for
prime sums on short multiplicative intervals.

## Importing arithmetic windows into the prime-power cone

Robert Polak · version 1.0 · 21 November 2026 · supplement to the 20 November preprint below

[PDF](papers/rh_theta_window_bridge.pdf) · [LaTeX source](papers/rh_theta_window_bridge.tex) · [Zenodo record (DOI 10.5281/zenodo.24415724)](https://zenodo.org/records/24415724)

How much control of the prime sum is enough to exclude a zeta zero? Fix
`0 ≤ s < 1/2`. If the integrated von Mangoldt residual satisfies
`F(x) ≤ C_s x^(1/2+s)` throughout `[Z, Z log Z]` for arbitrarily large `Z`,
with the same constant `C_s`, every nontrivial zero has real part at most
`1/2 + s`. The proof uses the positive kernel of the completed zeta function
to isolate a single zero, including its multiplicity.

The application to RH takes the arithmetic window estimate from Alpöge’s
work with Claude as an external input. The supplement supplies the analytic
step from the windows to zero exclusion and derives the determinant bound
`D_±(t) ≥ (1/32 − o(1)) t²` as `t → ∞` under RH.

The same residual gives an ordinary-prime test. Removing the higher prime powers from `F(4y) − 3F(2y) + 2F(y)` gives `P(y) − 2y log 2`, where `P(y)` is the dyadic tent sum over primes in `(y, 4y)`. For any fixed `δ > 0`, the bound `P(x/4) ≤ (x/2) log 2` throughout arbitrarily late windows `[X, X^(1+δ)]` implies RH. [bramble_trace’s criterion](https://rstation.ai/api/v1/posts/uw5cfb_J_Cj7XJO6/attachment?revision=24) and [Lee Crossley’s Lean formalisation](https://doi.org/10.5281/zenodo.24412819).

## An Explicit Prime-Power Cone Criterion for the Riemann Hypothesis

Robert Polak · version 1.2 · 20 November 2026

[Read the PDF](papers/Prime_Power_Cone_Criterion_en.pdf) · [Zenodo record and TeX source](https://doi.org/10.5281/zenodo.24406562)

The preprint builds on the August study by expressing the Riemann hypothesis
(RH) in terms of a determinant computed from two prime-power sums.

- **Criterion and bound.** RH holds exactly when the determinant immediately
  after each sufficiently large prime power is positive. Under RH, both
  one-sided versions are at least `0.006 t²` for `t = log x ≥ 11`. Interval
  arithmetic verifies positivity after all 6,148 prime powers up to 60,000.
- **Gaps between samples.** Let `p_m` be the first prime at or above `m⁴`.
  Divide the determinant by `√p_m`. For any fixed `0 < c ≤ 2B(0)` and block
  length `N`, the preprint gives an unconditional bound of
  `m exp(O_{c,N}((log log m)²))` on the starting index of the next block of
  `N` samples at or above `c`. For every fixed positive threshold, there are
  arbitrarily long blocks below it. RH is equivalent to a
  subpower bound on the wait for a positive sample. For every fixed `δ` in
  `(0, 1)`, a positive sample must occur within `H_δ m^δ` indices of every
  sufficiently large starting index `m`, for some constant `H_δ`.
- **Comparison with the Lean argument.** The paper identifies the shared
  von Mangoldt sum and explains what a theta-lift window estimate gives for
  the determinant. A bound for the signed combination of four shifted
  errors is still needed for the proposed short-interval transfer.

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
  Complete computational companion to the August study. It contains the finite
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

Run the non-destructive integrity audit for the August study from the
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
