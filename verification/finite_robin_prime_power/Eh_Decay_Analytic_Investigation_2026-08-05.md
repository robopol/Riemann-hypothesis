# Analytic investigation of the decay of the signed triangular discrepancy

Date: 2026-08-05

## 1. Bottom line

The downward cloud in the signed-triangular plot survives a full-population
finite audit; it is not created by the plotting sample. This investigation
produced three distinct results.

1. An exact unconditional prime-side reduction isolates the entire leading
   fluctuation:

   \[
    \boxed{
    E_h(q)=\frac{\Delta_H(q)}{2\sqrt q}
    +O_\lambda\!\left(\frac{(\log q)^3}{\sqrt q}\right),
    \qquad
    H=\lambda\sqrt q\log q.
    }
    \tag{1.1}
   \]

2. A full-population scan of all \(5{,}762{,}859\) prime-power events through
   \(10^8\) identifies the finite bulk scaling

   \[
    \boxed{
    \frac{|E_h(q)|}{Q(q)}
    \asymp_{\rm bulk}\frac{\log q}{\sqrt q}.
    }
    \tag{1.2}
   \]

3. Classical unconditional Selberg-integral estimates imply that the raw
   continuous numerator \(E^{\rm cont}_{h(x)}(\log x)\) tends to zero in
   normalized Lebesgue measure over dyadic blocks. This proves an almost-all
   continuous-center result, but it does not give a uniform bound at the
   prescribed sparse centers \(q=p^a\).

The remaining obstacle is no longer a hidden smooth integral. It is a
uniform, signed, state-adaptive short-interval estimate at every prime-power
center, together with the still-open central-secant/Turan branch.

## 2. Exact definition of the arithmetic fluctuation

Let \(q=p^a\), \(t=\log q\), and

\[
 h=\lambda\frac{\log q}{\sqrt q},
 \qquad
 H=qh=\lambda\sqrt q\log q.
\]

Define the odd additive triangular Chebyshev discrepancy

\[
 \boxed{
 \Delta_H(q)=
 \sum_{0<|r-q|<H}
 \operatorname{sgn}(r-q)
 \left(1-\frac{|r-q|}{H}\right)\Lambda(r).
 }
 \tag{2.1}
\]

Only prime powers contribute. The two halves are the triangularly weighted
von Mangoldt masses immediately to the right and left of \(q\). Their main
terms are the same, so their difference retains the cancellation discarded
by separate unsigned estimates.

The exact multiplicative version is

\[
\begin{aligned}
 \Delta_{\sinh}(q,h)=\frac2h\Bigg[&
 \sum_{q<r<qe^h}
 \Lambda(r)\sqrt{\frac qr}
 \sinh\!\left(\frac{h-\log(r/q)}2\right)\\
 &-\sum_{qe^{-h}<r<q}
 \Lambda(r)\sqrt{\frac qr}
 \sinh\!\left(\frac{h-\log(q/r)}2\right)
 \Bigg].
\end{aligned}
\tag{2.2}
\]

Using the half-weight central state of the regularized residual gives the
exact identity

\[
 \boxed{
 E_h(q)=
 \frac{\Delta_{\sinh}(q,h)}{2\sqrt q}
 +\kappa(h)\widehat b^\circ(q)
 +\sqrt q\,\beta(h)
 +\mathcal R_0(q,h),
 }
 \tag{2.3}
\]

where

\[
 \kappa(h)=\frac{2\sinh(h/2)}h-1,
 \qquad
 \beta(h)=\frac{2\sinh(h/2)}h-\cosh(h/2),
\]

and \(\mathcal R_0\) is the explicit trivial-zero remainder recorded in the
detailed derivation.

Elementary Chebyshev estimates give

\[
 \kappa(h)\widehat b^\circ(q)
 +\sqrt q\,\beta(h)+\mathcal R_0(q,h)
 =O_\lambda\!\left(\frac{(\log q)^2}{\sqrt q}\right).
 \tag{2.4}
\]

Expanding the hyperbolic kernel and replacing the multiplicative endpoints
by \(q\pm H\) gives

\[
 \Delta_{\sinh}(q,h)=\Delta_H(q)+O_\lambda((\log q)^3),
 \tag{2.5}
\]

which proves (1.1). Equations (2.3)--(2.5) are unconditional.

## 3. What elementary pointwise bounds prove

Brun--Titchmarsh for primes, together with an elementary count of higher
prime powers, gives

\[
 \sum_{qe^{-h}<r<qe^h}\Lambda(r)\ll H.
\]

Therefore

\[
 \boxed{
 |\Delta_H(q)|\ll H,
 \qquad
 |E_h(q)|\ll_\lambda\log q
 }
 \tag{3.1}
\]

uniformly and unconditionally. This is a valid theorem, but it loses a
factor of about \(q^{1/4}\) compared with the point cloud.

If an independent lower bound \(Q(q)\gg q^{1/4}\) held on a selected
subsequence, (3.1) would imply

\[
 \frac{|E_h(q)|}{Q(q)}
 \ll\frac{\log q}{q^{1/4}}\to0.
\]

Such a lower bound cannot be assumed uniformly. The exact radius satisfies

\[
 Q(q)^2=\mathfrak P(\log q)
 \left(2\sqrt q-\frac C2+\frac{\mathfrak P(\log q)}4+o(1)\right),
\]

so an anomalously small \(\mathfrak P\) also makes \(Q/q^{1/4}\) small.

## 4. Full-population scaling audit

For \(\lambda=1/4\), define

\[
 Z_1(q)=\frac{\sqrt q\,|E_h(q)|}{Q(q)\log q}.
\]

The audit recomputed all prime-power centers through \(10^8\), rather than
using only the \(17{,}836\)-point plotting sample.

| \(q\)-range | events | median \(Z_1\) | \(P_{90}\) | \(P_{99}\) | \(P_{99.9}\) | maximum |
|---|---:|---:|---:|---:|---:|---:|
| \(10^4\)--\(10^5\) | 8,420 | 0.2723 | 0.6717 | 1.0861 | 1.4541 | 1.8204 |
| \(10^5\)--\(10^6\) | 69,034 | 0.2951 | 0.7233 | 1.1312 | 1.4272 | 1.8489 |
| \(10^6\)--\(10^7\) | 586,400 | 0.2880 | 0.7097 | 1.1395 | 1.4888 | 1.9167 |
| \(10^7\)--\(10^8\) | 5,097,725 | 0.2933 | 0.7228 | 1.1412 | 1.4679 | 2.1662 |

Fitting

\[
 \frac{|E_h|}{Q}\asymp\frac{(\log q)^A}{\sqrt q}
\]

on logarithmic bins with \(q\ge10^5\) gives

\[
 A_{50}=0.968,
 \quad A_{90}=0.964,
 \quad A_{99}=0.986,
 \quad A_{99.9}=0.956.
\]

Thus \(A=1\) is a stable description of the bulk distribution. The maximum
of \(Z_1\) grows slowly. A descriptive four-decade fit adds approximately
\((\log q)^{0.316}\), so the observed maxima behave as though
\(A\approx1.32\). With only four maximum points, this is an extreme-value
diagnostic, not an asymptotic estimate. A proof target with \(A=3/2\) or
\(A=2\) is more realistic than demanding a uniform \(A=1\) bound.

The finite normalized collapse is not caused by sample selection:

- the raw plotting sample biases tail quantiles upward by about 11--24%;
- after the \(\sqrt q/\log q\) normalization, the bulk bias is mostly below
  5%;
- prime and higher-prime-power normalized distributions agree closely in
  the last decade;
- the sign split is 50.1568% positive and 49.8432% negative.

## 5. Why the exponent is natural

The window contains approximately

\[
 \frac{H}{\log q}\asymp\sqrt q
\]

primes. A prime contributes a von Mangoldt weight of order \(\log q\).
Square-root cancellation therefore predicts

\[
 |\Delta_H(q)|_{\rm typical}
 \asymp\sqrt{\sqrt q}\,\log q
 =q^{1/4}\log q.
\]

Substitution into (1.1) gives

\[
 |E_h(q)|_{\rm typical}\asymp\frac{\log q}{q^{1/4}}.
\]

When \(Q(q)\asymp q^{1/4}\), this becomes exactly

\[
 \frac{|E_h(q)|}{Q(q)}
 \asymp\frac{\log q}{\sqrt q}.
\]

The slow growth of the full-decade maxima is also natural: each decade
contains many more event centers, so a light-tailed cloud acquires an
extreme-value factor, heuristically of square-root-logarithmic size.

## 6. Spectral calibration under RH

Under RH, with the symmetric derivative convention,

\[
 E_h(t)=\sum_{\gamma>0}a_h(\gamma)\sin(\gamma t),
 \qquad
 a_h(\gamma)=
 \frac{2(\gamma-\sin(\gamma h)/h)}{\gamma^2+1/4}.
 \tag{6.1}
\]

This is a high-pass discrepancy filter. Fixed zeros with \(\gamma h\ll1\)
are suppressed quadratically, while for \(\gamma h\gg1\),
\(a_h(\gamma)\sim2/\gamma\). Thus the energy is generated primarily by the
moving band

\[
 \gamma\gtrsim1/h\asymp\frac{\sqrt q}{\log q},
\]

not by any fixed finite list of zeros.

Assuming RH and simple zeros, Besicovitch Parseval and the
Riemann--von Mangoldt formula give, for fixed \(h\to0\),

\[
 \boxed{
 \|E_h\|_{B^2}^2
 \sim\frac h6\log\frac1h.
 }
 \tag{6.2}
\]

The constant follows from

\[
 \int_0^\infty\frac{(y-\sin y)^2}{y^4}\,dy=\frac\pi6.
\]

Inserting \(h=\lambda(\log q)/\sqrt q\) identifies the natural scale

\[
 \boxed{
 |E_h|_{\rm natural}
 \asymp\sqrt{\frac\lambda{12}}
 \frac{\log q}{q^{1/4}}.
 }
 \tag{6.3}
\]

Equation (6.2) is rigorous for each fixed \(h\). The substitution in (6.3)
is a scale calibration, not by itself a theorem for a varying-window time
average and not a pointwise estimate.

About 98% of the limiting diagonal energy lies above \(\gamma=1/h\). This
explains why Taylor-expanding each fixed zero and then interchanging the
infinite sum cannot prove the eventwise law.

## 7. Unconditional almost-all decay

Put

\[
 h(x)=\lambda\frac{\log x}{\sqrt x},\qquad H(x)=xh(x).
\]

With the symmetric derivative at an atom, the continuous prime-side transfer
is

\[
 E^{\rm cont}_{h(x)}(\log x)
 =\frac1{2H(x)\sqrt x}\int_0^{H(x)}
 [\psi(x+u)+\psi(x-u)-2\psi^\circ(x)]\,du
 +O_\lambda\left(\frac{(\log x)^3}{\sqrt x}\right).
\]

After separating higher prime powers, define

\[
 \mathcal A_{H(x)}(x)=\frac1{2H(x)\sqrt x}\int_0^{H(x)}
 [\theta(x+u)+\theta(x-u)-2\theta^\circ(x)]\,du,
\]

where midpoint and right-continuous conventions differ only on a countable
set, and define

\[
 J(Y,u)=\int_Y^{2Y}|\theta(y+u)-\theta(y)-u|^2\,dy.
\]

For \(H=\lambda\sqrt X\log X\), Cauchy--Schwarz and Fubini give

\[
 \frac1X\int_X^{2X}|\mathcal A_{H(x)}(x)|^2\,dx
 \ll_\lambda\frac1{HX^2}\int_0^{C_\lambda H}
 [J(X/2,u)+J(X,u)]\,du.
 \tag{7.1}
\]

The Saffari--Vaughan Selberg-integral estimate is nontrivial uniformly from
\(u\ge X^{1/6+\varepsilon}\). Together with a trivial estimate below this
threshold, it yields, for every fixed \(0<\varepsilon<1/3\),

\[
 \boxed{
 \frac1X\int_X^{2X}
 |E^{\rm cont}_{h(x)}(\log x)|^2\,dx
 \ll_{\lambda,\varepsilon}
 (\log X)^2
 \exp\!\left[-c_\varepsilon
 \left(\frac{\log X}{\log\log X}\right)^{1/3}\right]
 +X^{-1+3\varepsilon}\log X
 +\frac{(\log X)^6}{X}.
 }
 \tag{7.2}
\]

The right side tends to zero. Therefore, for every \(\eta>0\),

\[
 \boxed{
 \frac1X\operatorname{meas}\{x\in[X,2X]:
 |E^{\rm cont}_{h(x)}(\log x)|>\eta\}\to0.
 }
 \tag{7.3}
\]

This is the strongest new analytic decay statement obtained in this
investigation. It still does not control the prime-power sequence. A
Lebesgue almost-all theorem allows a prescribed measure-zero set to be
entirely exceptional.

## 8. The exact open target

On a fixed positive state scale, the empirical law would follow from

\[
 \boxed{
 |\Delta_H(q)|\le C_A q^{1/4}(\log q)^A
 \quad\text{for every sufficiently large }q=p^a,
 }
 \tag{8.1}
\]

with any fixed \(A\). Because \(Q\) can be anomalously small, the robust
state-adaptive form must be stated for the exact discrepancy:

\[
 \boxed{
 |E_{h_\lambda(t)}(q)|
 \le C_{\lambda,A}Q(t)\frac{(\log q)^A}{\sqrt q},
 \qquad t=\log q.
 }
 \tag{8.2}
\]

The additive version \(|\Delta_H(q)|\ll Q(t)(\log q)^A\) also needs a
condition absorbing the \(O_\lambda((\log q)^3/\sqrt q)\) transfer error, or
one must work directly with the exact hyperbolic identity.

For the lower-cone induction, an absolute-value theorem is stronger than
necessary; only a one-sided lower-tail estimate is needed. It should ideally
be proved directly together with the secant term in

\[
 \mathfrak P(t+h)-\mathfrak P(t-h)
 +\mathcal T^\triangle(t,h)+2hQ(t)-hJ_q>0.
 \tag{8.3}
\]

The exact \(E_h\)-reduction does not prove the other open member
\(U_\lambda(q)\ge0\). The adaptive Turan inequality is itself an
RH-strength statement. Consequently, proving decay of \(E_h/Q\) would close
the discrepancy branch, but not the complete RH argument by itself.

At the present scale

\[
 T_0=1/h\asymp\frac{\sqrt q}{\log q},
 \qquad
 \frac{\log q}{\log T_0}=2+o(1).
\]

The analogous \(x=T^{2+o(1)}\) transition appears in the conjectural uniform
zero-exponential sums of Ford--Soundararajan--Zaharescu. Their results concern
different short-prime errors and do not imply (8.2). An average Selberg
estimate cannot simply be upgraded to a supremum over all event centers.

## 9. Best next analytic attacks

The most concrete routes are now:

1. prove a one-sided version of (8.1) with a safe logarithmic loss such as
   \(A=2\), first for prime centers and then for higher prime powers;
2. derive high-moment or large-deviation bounds for the signed triangular
   sum and a sampling theorem strong enough for the sparse event set;
3. attack the uniform zero exponential sum near \(q=T_0^{2+o(1)}\), where
   existing pair-correlation support stops short;
4. avoid a separate lower bound for \(Q\) by proving the combined
   state-correlated inequality (8.3) directly.

The fourth route is logically the cleanest for RH because it preserves the
correlation that separate numerator and denominator estimates lose near a
potentially small inherited \(\mathfrak P\)-state. Continuous-axis recurrence
does not by itself locate such a state at a prime-power event.

## 10. Reproducible artifacts

- [Exact prime-side derivation](Prime_Side_Eh_Decay_Derivation_2026-08-05.md)
- [Exact-decomposition regression test](test_prime_side_eh_decomposition.py)
- [Full scaling audit](Eh_Scaling_Law_Audit_2026-08-05.md)
- [Scaling audit script](audit_eh_scaling_law.py)
- [Full scaling report](eh_scaling_law_audit_report.json)
- [Scaling smoke tests](test_audit_eh_scaling_law.py)
- [Spectral and Selberg analysis](Spectral_Eh_Decay_Analysis_2026-08-05.md)

Run the fast regression suite from the repository root:

```powershell
python -m unittest test_prime_side_eh_decomposition.py -v
python -m unittest test_audit_eh_scaling_law.py -v
```

Both suites pass.  The exact reduction, scaling audit, and rigor boundary are
also incorporated into the standalone findings paper
`papers/Finite_Robin_Verification_via_CA_Prime_Power_Reduction_en.tex`.

## 11. Primary literature boundary

- B. Saffari and R. C. Vaughan, Section 6,
  [On the fractional parts of x/n and related sequences. II](https://www.numdam.org/item/AIF_1977__27_2_1_0/).
- H. L. Montgomery and K. Soundararajan,
  [Primes in short intervals](https://arxiv.org/abs/math/0409258).
- T. H. Chan,
  [More precise pair correlation of zeros and primes in short intervals](https://arxiv.org/abs/math/0206292).
- K. Ford, K. Soundararajan, and A. Zaharescu,
  [On the distribution of imaginary parts of zeros of the Riemann zeta function, II](https://arxiv.org/abs/0805.2745).

These sources support the averaged/pair-correlation boundary. None supplies
the state-adaptive, every-prime-power inequality (8.2) or (8.3).
