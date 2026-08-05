# Finite scaling audit for the signed triangular discrepancy

Date: 2026-08-05

## 1. Question and conclusion

This audit tests the finite-range law suggested by the upper-right panel of
the signed-triangular figure.  At every prime-power event \(q=p^a\), put

\[
 h=\frac{\log q}{4\sqrt q},
 \qquad
 E_h=\bar p-\frac{\Phi(t+h)-\Phi(t-h)}{2h},
 \qquad t=\log q,
\]

and let \(Q\) denote the cone radius called \(R\) in the derivation notes.
The tested normalization is

\[
 \boxed{
 Z_A(q)=\frac{\sqrt q\,|E_h(q)|}
              {Q(q)(\log q)^A}.
 }
\]

The main finite finding is stronger than the visual impression from the
17,836-point plotting sample:

\[
 \boxed{
 \frac{|E_h(q)|}{Q(q)}
 \text{ has a bulk scale very close to }
 \frac{\log q}{\sqrt q}
 \quad(10^4\le q\le10^8).
 }
\]

This conclusion was recomputed on **all 5,762,859 prime-power events**, not
on the plotting sample.  For \(A=1\), the median, 90th, 99th, and 99.9th
percentiles of \(Z_1\) are nearly stationary across the last four complete
decades.  Their across-decade coefficients of variation are only 3.12%,
2.98%, 2.00%, and 1.53%, respectively.

The maxima require a separate qualification.  The full-decade maximum of
\(Z_1\) rises mildly from 1.820 to 2.166.  A four-point descriptive fit is
compatible with an additional very slow extreme-value factor, approximately
\((\log q)^{0.316}\).  Consequently the finite data suggest

- \(A\simeq1\) for the bulk distribution;
- an effective \(A\simeq1.3\) for the observed maxima;
- a target such as \(A=3/2\) or \(A=2\) would be a conservative analytic
  envelope to try.

None of these finite fits is a proof of a universal bound.  In particular,
the computation supplies no data for \(q>10^8\).

## 2. Population, computation, and fit protocol

The companion script is
`audit_eh_scaling_law.py`.  It independently performs the aggregation and
sampling audit, while importing the already audited prime-power support and
state evaluators from the signed-triangular scanner.  Run it from the
repository root with

```powershell
python audit_eh_scaling_law.py `
  --limit 100000000 --lambda-value 0.25 `
  --segment-span 10000000 --chunk-size 250000 `
  --sample-count 30000 --report eh_scaling_law_audit_report.json
```

It writes `eh_scaling_law_audit_report.json`.  The full calculation took
54.42 seconds in the recorded Python 3.13 / NumPy 2.3.5 environment.

The event population is:

| population | count |
|---|---:|
| primes \(q=p\) | 5,761,455 |
| higher prime powers \(q=p^a,\ a\ge2\) | 1,404 |
| all prime powers | 5,762,859 |
| \(E_h>0\) | 2,890,468 |
| \(E_h<0\) | 2,872,391 |
| binary64 zeros | 0 |

To estimate \(A\), every \(q\)-decade is divided into ten equal-width bins
in \(\log q\).  In each bin the relevant quantile of

\[
 \sqrt q\,\frac{|E_h|}{Q}
\]

is computed.  Its logarithm is then regressed on the logarithm of
\(\log q\).  Under the proposed model the fitted slope is \(A\).  Equal log
bins prevent the last, densest decade from dominating the fit merely because
it contains more primes.

The reported standard errors are ordinary regression diagnostics.  They do
not model arithmetic dependence between neighboring prime-power windows and
must not be read as rigorous confidence intervals.

## 3. Full-population decadal quantiles

The following table is based on all events in each decade.  It also shows
that the fall is not produced only by a growing denominator: the median of
\(|E_h|\) itself decreases throughout the tail.

| \(q\)-range | events | median \(|E_h|\) | median \(|E_h|/Q\) | median \(Z_1\) | \(P_{90}(Z_1)\) | \(P_{99}(Z_1)\) | \(P_{99.9}(Z_1)\) | max \(Z_1\) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| \(10^4\)--\(10^5\) | 8,420 | 0.060920 | 0.013541 | 0.272313 | 0.671667 | 1.086138 | 1.454088 | 1.820441 |
| \(10^5\)--\(10^6\) | 69,034 | 0.044177 | 0.005568 | 0.295081 | 0.723275 | 1.131155 | 1.427215 | 1.848881 |
| \(10^6\)--\(10^7\) | 586,400 | 0.028322 | 0.002003 | 0.287972 | 0.709653 | 1.139468 | 1.488830 | 1.916707 |
| \(10^7\)--\(10^8\) | 5,097,725 | 0.018836 | 0.000746 | 0.293324 | 0.722770 | 1.141239 | 1.467949 | 2.166157 |

Thus the normalized bulk is exceptionally stable.  Across these four rows,
the relative ranges (maximum minus minimum, divided by the mean) are 7.93%
for the median, 7.30% for \(P_{90}\), 4.90% for \(P_{99}\), and 4.22% for
\(P_{99.9}\).

## 4. Estimate of the logarithmic exponent

The fitted exponents are:

| fitted population | quantile | \(A\), \(q\ge10^4\) | \(A\), \(q\ge10^5\) |
|---|---:|---:|---:|
| all events | 50% | \(1.130\pm0.060\) | \(0.968\pm0.095\) |
| all events | 90% | \(1.129\pm0.060\) | \(0.964\pm0.087\) |
| all events | 99% | \(1.186\pm0.070\) | \(0.986\pm0.094\) |
| all events | 99.9% | \(1.286\pm0.085\) | \(0.956\pm0.115\) |

The \(q\ge10^4\) fit still contains a visible transition from the earlier
range.  Once the cutoff is moved to \(10^5\), all four quantile slopes cluster
between 0.956 and 0.986.  Together with the stationary decadal table, this
makes

\[
 \boxed{A=1}
\]

the best parsimonious finite-range description of the bulk.  The data do not
justify treating 1.130 as a stable asymptotic exponent.

## 5. Signs and higher prime powers

The global sign split is 50.1568% positive and 49.8432% negative.  In the
last four decades the positive fractions are 50.095%, 50.481%, 50.216%, and
50.146%.  The fall therefore does not come from an eventual fixed sign.

At \(10^7\le q<10^8\), the normalized distributions separated by type and
sign are:

| group | events | median \(Z_1\) | \(P_{90}\) | \(P_{99}\) | \(P_{99.9}\) | max |
|---|---:|---:|---:|---:|---:|---:|
| primes | 5,096,876 | 0.293324 | 0.722771 | 1.141249 | 1.467930 | 2.166157 |
| higher prime powers | 849 | 0.296766 | 0.718863 | 1.121019 | 1.618901 | 1.696293 |
| \(E_h>0\) | 2,556,303 | 0.292343 | 0.720275 | 1.131964 | 1.448687 | 2.166157 |
| \(E_h<0\) | 2,541,422 | 0.294310 | 0.725356 | 1.150190 | 1.484951 | 1.970255 |

The prime and higher-power bulk distributions agree surprisingly well in
the last decade.  The higher-power population is nevertheless too small and
heterogeneous for a stable exponent fit: only 1,404 such events occur in the
entire scan, and their fitted \(R^2\) values are much weaker.  No separate
higher-power asymptotic law should be claimed from this dataset.

The full \(Z_1\) maxima in the last four decades occur at

| decade | max \(Z_1\) | witness \(q\) | type | sign of \(E_h\) |
|---|---:|---:|---|---|
| \(10^4\)--\(10^5\) | 1.820441 | 63,863 | prime | negative |
| \(10^5\)--\(10^6\) | 1.848881 | 250,673 | prime | positive |
| \(10^6\)--\(10^7\) | 1.916707 | 5,579,347 | prime | positive |
| \(10^7\)--\(10^8\) | 2.166157 | 34,763,537 | prime | positive |

Thus the rescaled extremes are not being driven by higher prime powers in
these four decades.

## 6. Extreme-value diagnostic

For \(A=1\), a regression using only the four displayed decade maxima gives

\[
 \max Z_1 \asymp (\log q)^{0.316\pm0.119},
 \qquad R^2=0.778.
\]

Equivalently, regressing the same maxima on the number \(N\) of events in a
decade gives

\[
 \max Z_1 \asymp N^{0.0262\pm0.0081},
 \qquad R^2=0.839.
\]

Only four points are present, so neither expression is a reliable
extrapolation.  They do show why the bulk exponent and a deterministic
uniform exponent should not be conflated.  The 99.9th percentile remains
between 1.427 and 1.489, while the maximum rises as each new decade supplies
roughly ten times as many opportunities for a rare event.  This behavior is
compatible with a light-tailed extreme-value correction such as a power of
\(\log q\), but the data cannot determine its ultimate form.

The finite maximum fit corresponds to the descriptive envelope

\[
 \frac{|E_h|}{Q}
 \lesssim \frac{(\log q)^{1.32}}{\sqrt q}
\]

over the checked tail.  An attempted proof should use a rounder and safer
target, for example

\[
 \boxed{
 |E_h(q)|\le C Q(q)\frac{(\log q)^{3/2}}{\sqrt q}
 }
\]

or even exponent \(2\).  Any fixed logarithmic exponent would still imply
\(|E_h|/Q\to0\).  The numerical fit does not supply the required universal
constant \(C\).

## 7. Plot-sample selection bias

The plotted blue/gold cloud is deliberately not a random sample.  The
generator includes

1. every event through \(10^4\);
2. every higher prime power through \(10^8\);
3. event indices nearest to 30,000 geometrically spaced targets.

Duplicate targets collapse, leaving 17,836 plotted events.  The resulting
tail sampling fractions are:

| decade | full events | plotted events | sampled fraction | sample bias in median \(Z_1\) | sample bias in \(P_{99}(Z_1)\) | fraction of full max captured |
|---|---:|---:|---:|---:|---:|---:|
| \(10^4\)--\(10^5\) | 8,420 | 3,578 | 42.494% | -0.46% | -0.56% | 100.00% |
| \(10^5\)--\(10^6\) | 69,034 | 4,018 | 5.820% | -1.05% | +5.03% | 90.35% |
| \(10^6\)--\(10^7\) | 586,400 | 4,214 | 0.719% | +0.48% | +1.88% | 99.92% |
| \(10^7\)--\(10^8\) | 5,097,725 | 4,746 | 0.093% | -2.71% | -0.76% | 90.03% |

The raw \(|E_h|/Q\) sample quantiles are biased upward by roughly 11--24% in
the last four decades because log-uniform sampling gives relatively more
weight to the lower-\(q\) part of each decade.  After the proposed
\(\sqrt q/\log q\) normalization, the central and 99th-percentile biases are
mostly within 5%.  Therefore the observed distributional collapse is not an
artifact of the plotting sample.

All 1,404 higher prime powers are included in the figure.  They form 7.87%
of plotted points but only 0.0244% of the full population, an overrepresentation
factor of about 323.  The gold rings are useful for visibility but must not be
read as their population frequency.  The orange line in the figure uses the
archived full-scan decade maxima, so it is not a sample maximum.

## 8. Validation and limitations

The following consistency checks passed:

- the recomputed prime, higher-power, and total event counts equal the
  independently audited population counts;
- every recomputed decade maximum of \(|E_h|/Q\) agrees exactly in binary64
  with the archived \(10^8\) report;
- the reconstructed 17,836 sample has the same event count and the same
  quantile summaries as the exported plotting CSV;
- no binary64 zero of \(E_h\) occurs.

This script is an independent aggregation and selection-bias audit, not an
independent derivation of the prime-side state.  Formula and precision QA are
implemented by `test_exploratory_signed_triangular_direct_qa.py`, including an
80-digit direct calculation with mpmath.  That audit found a largest
high-precision spot-check change in \(E_h\) of about
\(6.83\times10^{-7}\), whereas the full last-decade median \(|E_h|\) here is
0.01884.  The bulk scaling is therefore far above the observed numerical
error scale, although fitted exponents remain non-certified binary64
statistics.

The decisive limitation is logical rather than numerical: stable quantiles,
or even one hundred million verified events, do not bound the supremum in
all later decades.  A proof still needs a deterministic prime-side estimate
valid for every sufficiently large prime power.

## 9. Reproducibility hashes

- `audit_eh_scaling_law.py` SHA-256:
  `f8c3d03963705cddead22ab07ff67a8f8d17855ee2d363cb8d0f5d0eebe6100b`;
- `eh_scaling_law_audit_report.json` SHA-256:
  `0aaf4eb5dd456e94417853dcf4609d43c0113fb5b14a097928589e9f56d3fa4b`;
- `test_audit_eh_scaling_law.py` SHA-256:
  `f2987d0407dd59af4d0606d12f617ba50ff7e0afb9f587d78c12116226e35d6d`;
- archived full scan SHA-256:
  `b464e154ee46a8caf9a223dda8c60ebcacd17eefa4d9548c33b10e16232743a4`;
- plotted sample CSV SHA-256:
  `e35fd4e97c829f5310bdba29d5d668a83f9bcb7b1d7ed40e9a2c695b9bd218bb`.

## 10. Bottom line

The downward motion of the cloud is real and survives a full-population
audit.  The strongest evidence-supported finite statement is

\[
 \boxed{
 Z_1(q)=\frac{\sqrt q\,|E_h(q)|}{Q(q)\log q}
 \text{ has an approximately stationary bulk distribution through }10^8.
 }
\]

What remains open is upgrading that finite distributional law to a
deterministic all-event envelope.  The empirical target is now much sharper:
prove a bound with the square-root factor and any fixed polylogarithmic loss,
rather than trying to infer an unspecified decay directly from the plot.
