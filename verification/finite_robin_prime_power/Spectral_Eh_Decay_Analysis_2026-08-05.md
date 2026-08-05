# Spectral analysis of the signed triangular discrepancy

Date: 2026-08-05

## 1. Scope and conclusion

This note analyzes the finite scaling observed through \(q\le10^8\):

\[
 \frac{|E_h(q)|}{Q(q)}\approx_{\rm bulk} \frac{\log q}{\sqrt q},
 \qquad
 h=\lambda\frac{\log q}{\sqrt q},
\]

where

\[
 E_h(t)=\bar p(t)-
 \frac{\mathfrak P(t+h)-\mathfrak P(t-h)}{2h}
 =\frac{\mathcal T^\triangle(t,h)}{2h}.
\]

There are four logically different statements.

1. **Unconditional pointwise bound.** For fixed \(\lambda>0\), elementary
   short-interval estimates give
   \[
    E_{h(t)}(t)=O_\lambda(t),\qquad t=\log x.
   \]
   This does not show decay.
2. **Unconditional almost-all decay.** The Saffari--Vaughan Selberg-integral
   estimate at \(H=xh\asymp\sqrt x\log x\) implies
   \[
    E^{\rm cont}_{h(x)}(\log x)\longrightarrow0
   \]
   in normalized Lebesgue measure on dyadic \(x\)-blocks. This is an analytic
   continuous-center counterpart to the observed fall, but it says nothing at
   the prescribed measure-zero set \(x=p^a\).
3. **RH-conditional natural scale.** For each fixed \(h\), the zero-side
   triangular filter has
   diagonal spectral energy
   \[
    \|E_h\|_{B^2}^2\sim \frac h6\log\frac1h
   \]
   if RH and simplicity of the zeros are assumed. Substitution of the
   adaptive \(h\) only calibrates the predicted scale
   \[
    |E_h|_{\rm predicted}\asymp
    \sqrt{\frac{\lambda}{12}}\,
    \frac{\log x}{x^{1/4}}.
   \]
   It is not a varying-path RMS theorem. When \(\mathfrak P(t)\) is on a
   fixed positive scale, the corresponding predicted ratio has scale
   \((\log x)/\sqrt x\).
4. **Uniform every-event target.** None of the preceding statements proves
   this scaling at every \(q=p^a\). Such a statement requires a new uniform
   signed short-interval estimate, or a sufficiently strong route through
   zero-exponential sums, in the
   \(\alpha=\log x/\log(1/h)=2+o(1)\) regime. This is on the conjectural side
   of strong pair correlation.

Thus the finite downward cloud survives the full-population audit. Its scale
is consistent with the fixed-window spectral calibration, and numerator
decay is proved for almost all continuous centers. The missing step is a
uniform theorem at every prime-power event, coupled to the adaptive
denominator and the central-secant branch.

## 2. Exact RH-side filter

Assume RH. The absolutely convergent positive expansion is

\[
 \mathfrak P(t)=2\sum_{\gamma>0}
 \frac{1-\cos(\gamma t)}{\gamma^2+1/4}.
\]

At a prime-power event, the derivative is interpreted by the symmetric
explicit-formula convention. Formally differentiating in that convention and
subtracting the central secant gives

\[
 \boxed{
 E_h(t)=\sum_{\gamma>0}a_h(\gamma)\sin(\gamma t),
 \qquad
 a_h(\gamma)=
 \frac{2\left(\gamma-\sin(\gamma h)/h\right)}{\gamma^2+1/4}.
 }
\]

The series is not absolutely convergent: \(a_h(\gamma)\sim2/\gamma\) at
high frequency. It is square summable and therefore defines a Besicovitch
\(B^2\) trigonometric series. This distinction is essential; ordinary
termwise absolute summation would lose the cancellation under investigation.

The response of a general zero \(\rho\), with
\(\mu=\rho-1/2\), is

\[
 E_{\rho,h}(t)=
 \frac{e^{\mu t}}{\rho(1-\rho)}
 \left(\frac{\sinh(\mu h)}h-\mu\right).
\]

For a fixed zero and \(|\mu|h\ll1\),

\[
 E_{\rho,h}(t)=
 \frac{\mu^3h^2}{6\rho(1-\rho)}e^{\mu t}
 +O_\rho(h^4e^{\Re(\mu)t}).
\]

Consequently, every fixed off-critical zero is strongly suppressed by the
shrinking window. This is the spectral form of the no-go observation: decay
of \(E_h\), even absolute decay, does not by itself exclude a growing
off-line mode.

## 3. Low/high-frequency split

Put \(Y=1/h\) and \(y=\gamma h\). Then

\[
 a_h(\gamma)=
 2h\frac{y-\sin y}{y^2+h^2/4}.
\]

For \(\gamma h\ll1\),

\[
 a_h(\gamma)
 =\frac{\gamma^3}{3(\gamma^2+1/4)}h^2
 +O\left(\frac{\gamma^5}{\gamma^2+1/4}h^4\right),
\]

so each fixed low mode is killed quadratically. For \(\gamma h\gg1\),

\[
 a_h(\gamma)=\frac2\gamma
 +O\left(\frac1{h\gamma^2}\right).
\]

Thus the filter is a high-pass discrepancy filter, not a smoothing filter.
The observed size is generated mainly by zeros whose height moves upward as

\[
 \gamma\gtrsim\frac1h
 \asymp\frac{\sqrt x}{\log x}.
\]

This explains why a Taylor expansion for each fixed zero cannot prove a
uniform result. The expansion is nonuniform in the infinite zero sum.

## 4. Exact diagonal energy and the predicted power law

Assume RH and simple zeros. Besicovitch Parseval gives, for fixed \(h\),

\[
 \lim_{L\to\infty}\frac1L\int_0^L |E_h(t)|^2\,dt
 =\frac12\sum_{\gamma>0}a_h(\gamma)^2.
\]

Using the Riemann--von Mangoldt formula as \(h\downarrow0\),

\[
 \sum_{\gamma>0}a_h(\gamma)^2
 \sim\frac{2h}{\pi}\log\frac1h
 \int_0^\infty\left(\frac{y-\sin y}{y^2}\right)^2dy.
\]

The remaining integral is exact:

\[
 \int_0^\infty\frac{(y-\sin y)^2}{y^4}\,dy=\frac\pi6.
\]

Indeed,

\[
 \frac{y-\sin y}{y^2}
 =\int_0^1\frac{1-\cos(uy)}y\,du,
\]

and

\[
 \int_0^\infty
 \frac{(1-\cos uy)(1-\cos vy)}{y^2}\,dy
 =\frac\pi2\min(u,v).
\]

Therefore

\[
 \boxed{
 \sum_{\gamma>0}a_h(\gamma)^2
 \sim\frac h3\log\frac1h,
 \qquad
 \|E_h\|_{B^2}^2
 \sim\frac h6\log\frac1h.
 }
\]

If multiple zeros are allowed, equal ordinates must be combined before
Parseval. RH alone gives the same lower diagonal scale and an upper bound with
at most an extra logarithm using the classical multiplicity bound. The clean
constant above uses simplicity.

Only as a two-limit calibration, first taking the long-time mean for fixed
\(h\) and then sending \(h\downarrow0\), substitute

\[
 h=\lambda t e^{-t/2},\qquad x=e^t.
\]

Since \(\log(1/h)=t/2+O(\log t)\),

\[
 \boxed{
 |E_h|_{\rm diagonal\ scale}
 \asymp\sqrt{\frac{\lambda}{12}}\,
 t e^{-t/4}
 =\sqrt{\frac{\lambda}{12}}\,
 \frac{\log x}{x^{1/4}}.
 }
\]

In the RH regime,

\[
 Q(t)=\sqrt{2e^{t/2}\mathfrak P(t)}\,(1+o(1)).
\]

Whenever \(\mathfrak P(t)\) remains on a fixed positive scale, division gives

\[
 \boxed{
 \frac{|E_h|_{\rm diagonal\ scale}}{Q(t)}
 \asymp
 \sqrt{\frac{\lambda}{24\mathfrak P(t)}}
 \frac{\log x}{\sqrt x}.
 }
\]

This fixed-window calculation calibrates the stationary rescaling seen in
the finite data.  It is not a norm or theorem along the varying path
\(h(t)=\lambda t e^{-t/2}\), not a local or dyadic RMS theorem, and not a
pointwise envelope.

## 5. Numerical check of the spectral energy

This optional spot-check uses an external Odlyzko zero table from the broader
working archive.  The table is not part of the compact public companion
package, and none of the analytic statements below depends on it.

The first 100,000 stored zero ordinates, ending at
\(\gamma=74920.827498994\), give the following partial diagonal energies for
\(\lambda=1/4\):

| \(x\) | \(h\) | \(1/h\) | \(\sum_{\gamma\le74920}a_h(\gamma)^2/[h\log(1/h)/3]\) | fraction from \(\gamma h\le1\) |
|---:|---:|---:|---:|---:|
| \(10^4\) | 0.023026 | 43.43 | 0.9451 | 0.0088 |
| \(10^5\) | 0.009102 | 109.87 | 0.9522 | 0.0091 |
| \(10^6\) | 0.003454 | 289.53 | 0.9520 | 0.0107 |
| \(10^7\) | 0.001274 | 784.78 | 0.9395 | 0.0119 |
| \(10^8\) | 0.000461 | 2171.47 | 0.8996 | 0.0133 |

The falling ratio at \(10^8\) is the omitted tail above the last stored zero.
The theoretical asymptotic fraction below \(\gamma h=1\) is only

\[
 \frac{\displaystyle\int_0^1((y-\sin y)/y^2)^2dy}{\pi/6}
 =0.016659\ldots.
\]

Thus roughly 98% of the diagonal energy is produced above the moving cutoff
\(1/h\). This directly confirms the nonuniform-summation obstruction.

## 6. Prime-side form and an unconditional pointwise bound

Let \(H=xh=\lambda\sqrt x\log x\). The atomic part of the exact triangular
identity satisfies

\[
 \frac{|\mathcal T^\triangle_{\rm atom}(t,h)|}{2h}
 \le\frac12\sum_{xe^{-h}<r<xe^h}\frac{\Lambda(r)}{\sqrt r}.
\]

Brun--Titchmarsh for the prime contribution, plus an elementary count of
higher prime powers, gives

\[
 \sum_{xe^{-h}<r<xe^h}\Lambda(r)\ll H,
\]

and hence the atomic term is \(O_\lambda(\log x)\). Chebyshev bounds in the
exact prime-side formula give

\[
 |\mathfrak P'(t)|+|\mathcal F'(t)|\ll\sqrt x\log x.
\]

The odd smooth contribution is therefore

\[
 \frac1{2h}\left|
 \int_0^h(h-s)[\mathcal F(t+s)-\mathcal F(t-s)]\,ds
 \right|
 \ll h^2\sqrt x\log x
 \ll\frac{(\log x)^3}{\sqrt x}.
\]

Consequently,

\[
 \boxed{E_{h(t)}(t)=O_\lambda(\log x)}
\]

unconditionally and uniformly. This is far weaker than the observed
\((\log x)x^{-1/4}\) size because it takes absolute values across the two
sides of the window.

## 7. Unconditional almost-all numerator decay

There is nevertheless an unconditional averaged result at exactly the target
window length.  Define the continuous-center version by

\[
 h(x)=\lambda\frac{\log x}{\sqrt x},\qquad
 H(x)=xh(x),
\]

and use the symmetric derivative at an atom.  The continuous prime-side
transfer is

\[
 E^{\rm cont}_{h(x)}(\log x)
 =\frac1{2H(x)\sqrt x}\int_0^{H(x)}
 [\psi(x+u)+\psi(x-u)-2\psi^\circ(x)]\,du
 +O_\lambda\left(\frac{(\log x)^3}{\sqrt x}\right).
\]

After separating higher prime powers, the prime contribution uses
\(\theta^\circ(x)\), the midpoint value at a prime.  Replacing it by the
standard right-continuous value affects only a countable set.  Put

\[
 \mathcal A_{H(x)}(x)=\frac1{2H(x)\sqrt x}\int_0^{H(x)}
 [\theta(x+u)+\theta(x-u)-2\theta^\circ(x)]\,du
\]

and

\[
 J(Y,u)=\int_Y^{2Y}
 |\theta(y+u)-\theta(y)-u|^2\,dy.
\]

For \(H=\lambda\sqrt X\log X\), Cauchy--Schwarz, Fubini, and
\(H(x)\asymp_\lambda H\) on \([X,2X]\) give

\[
 \frac1X\int_X^{2X}|\mathcal A_{H(x)}(x)|^2dx
 \ll_\lambda\frac1{HX^2}\int_0^{C_\lambda H}
 [J(X/2,u)+J(X,u)]\,du.
\]

Saffari and Vaughan, Section 6, prove that for every fixed
\(0<\epsilon<1/3\), uniformly for
\(X^{1/6+\epsilon}\le u\le C_\lambda H\),

\[
 J(X/2,u)+J(X,u)\ll_\epsilon
 u^2X\exp\left[-c_\epsilon
 \left(\frac{\log X}{\log\log X}\right)^{1/3}\right].
\]

The trivial estimate below \(X^{1/6+\epsilon}\) and the transfer remainder
then yield

\[
 \boxed{
 \frac1X\int_X^{2X}
 |E^{\rm cont}_{h(x)}(\log x)|^2dx
 \ll_{\lambda,\epsilon}
 (\log X)^2
 \exp\left[-c_\epsilon
 \left(\frac{\log X}{\log\log X}\right)^{1/3}\right]
 +X^{-1+3\epsilon}\log X
 +\frac{(\log X)^6}{X}.
 }
\]

The right side tends to zero.  Chebyshev's inequality therefore gives, for
every \(\eta>0\),

\[
 \boxed{
 \frac1X\operatorname{meas}\{x\in[X,2X]:
 |E^{\rm cont}_{h(x)}(\log x)|>\eta\}\to0.}
\]

This is the strongest unconditional decay obtained here. It does **not**
control the values at prime powers. A Lebesgue \(L^2\) theorem permits a
measure-zero prescribed sequence to lie entirely in its exceptional set.

Primary source:

- B. Saffari and R. C. Vaughan, Section 6, *On the fractional parts of
  \(x/n\) and related sequences. II*, Annales de l'Institut Fourier 27
  (1977), 1--30,
  [DOI 10.5802/aif.649](https://doi.org/10.5802/aif.649).

## 8. Why the fixed-window calibration and every-event maximum remain open

The sharp short-interval variance model is

\[
 \int_X^{2X}|\psi(x+H)-\psi(x)-H|^2dx
 \sim XH\log\frac XH.
\]

At \(H=\lambda\sqrt X\log X\), this model would predict

\[
 \operatorname{RMS}(E_h)\asymp
 \frac{\log X}{X^{1/4}},
\]

matching the fixed-window spectral calibration. Under RH, this variance is
linked to strong pair correlation. T. H. Chan's Theorem 2.1 gives a precise
quantitative version of the relation between the pair-correlation function
and the short-interval second moment. At the present cutoff

\[
 T_0=\frac1h\asymp\frac{\sqrt X}{\log X},
 \qquad
 \frac{\log X}{\log T_0}=2+o(1).
\]

This is the \(\alpha\simeq2\) side, beyond the classically proved
\(\alpha\le1\) support.

Primary source:

- T. H. Chan, Theorem 2.1, *More precise pair correlation of zeros and
  primes in short intervals*, J. London Math. Soc. 68 (2003), 579--598,
  [arXiv:math/0206292](https://arxiv.org/abs/math/0206292).

Ford--Soundararajan--Zaharescu isolate the still stronger uniform obstruction.
Their Conjecture 3 asks, under RH, for

\[
 \sum_{0<\gamma\le T}x^{i\gamma}=o(T)
\]

uniformly from \(T^2/(\log T)^5\le x\le T^A\). Our relation
\(x=T_0^{2+o(1)}\) lies in the analogous transition. Their Theorem 1.3 says
that Conjecture 3 implies their Conjectures 1 and 2 about fractional-part
distributions; it does not prove our signed triangular estimate. A later
theorem in the same paper relates a different short-prime conjecture to
Conjecture 3. In the reverse direction, Conjecture 3 yields only a
\(o(\sqrt x\log x)\) error for the relevant prime increment, much coarser
than the \(x^{1/4}\) scale required here. The comparison is therefore a
spectral-regime analogy, not an implication for our kernel.

Primary source:

- K. Ford, K. Soundararajan and A. Zaharescu, Conjecture 3 and the two
  following theorems,
  *On the distribution of imaginary parts of zeros of the Riemann zeta
  function, II*, Math. Ann. 343 (2009), 487--505,
  [arXiv:0805.2745](https://arxiv.org/abs/0805.2745).

Even the sharp pair-correlation variance is an average, not a supremum. To
prove the plotted full-decade maximum continues to fall, one needs one of:

1. a uniform exponential-sum bound of the Ford--Soundararajan--Zaharescu
   type with enough quantitative saving;
2. high-moment/large-deviation estimates plus a rigorous sampling theorem for
   all prime-power centers;
3. a direct prime-side signed triangular inequality that exploits the event
   correlation and bypasses separate numerator/denominator estimates.

## 9. The denominator obstruction

The graph plots \(|E_h|/Q\), not merely \(|E_h|\). Under RH,

\[
 \mathfrak P(t)=2\sum_{\gamma>0}
 \frac{1-\cos(\gamma t)}{\gamma^2+1/4}\ge0,
 \qquad
 Q^2=\frac{\mathfrak P^2}{4}+2k\mathfrak P.
\]

The fixed-window diagonal calibration assumes that \(\mathfrak P(t)\) is on
a normal positive scale. It supplies no uniform lower bound for
\(\mathfrak P(\log q)\). Under RH, almost-periodic recurrence gives small
values on the continuous \(t\)-axis, but it does not place them at
\(t=\log(p^a)\). An eventwise lower bound is itself a nontrivial
Diophantine/zero-distribution problem. Therefore even a uniform theorem
\(E_h=o(1)\) would not by itself prove the displayed normalized envelope.

The proof target remains the combined signed inequality

\[
 \mathfrak P(t+h)-\mathfrak P(t-h)
 +\mathcal T^\triangle(t,h)+2hQ(t)-hJ_q>0,
\]

or the two correlated lemmas \(U_\lambda(q)\ge0\) and
\(M_\lambda(q)>0\). Separating \(E_h\) and \(Q\) loses precisely the
correlation needed near a small-\(\mathfrak P\) recurrence.

## 10. Strongest justified takeaway

The visible downward motion has both an analytic scale calibration and an
unconditional continuous-center mean-square counterpart:

\[
 \boxed{
 \text{fixed-}h\text{ diagonal energy predicts the scale}
 \quad
 \frac{\log x}{x^{1/4}}.}
\]

Moreover, classical unconditional Selberg-integral technology proves that the
raw numerator tends to zero for almost all continuous centers at this window
length. What it does not prove is the statement required by the RH mechanism:

\[
 \boxed{
 |E_{h_\lambda(t)}(p^a)|
 \le C_{\lambda,A}\,Q(t)\frac{(\log(p^a))^A}{\sqrt{p^a}},
 \qquad t=\log(p^a),
 \quad\text{for every sufficiently large prime power }p^a,}
\]

for some fixed exponent \(A\) (the bulk suggests \(A=1\), while a uniform
extreme-value envelope may need a larger logarithmic power). The empirical
law is observed and consistent with the fixed-window calibration, but its
uniform eventwise version is a strong signed short-interval problem coupled
to the state-dependent denominator, not a consequence of the finite graph.
