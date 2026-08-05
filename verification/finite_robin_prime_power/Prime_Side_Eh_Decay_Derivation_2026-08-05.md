# Prime-side derivation of the observed decay of the signed triangular discrepancy

**Date:** 2026-08-05  
**Scope:** derive an exact prime-side formula for \(E_h\), isolate the term
responsible for the observed decay, and state precisely what remains to be
proved. This note does not assume RH unless explicitly stated.

## 1. Result in one line

Let \(q=p^a\), \(t=\log q\), and

\[
 h=\lambda\frac{\log q}{\sqrt q},\qquad \lambda>0\ \text{fixed}.
\]

Define \(H=qh=\lambda\sqrt q\log q\) and the additive odd triangular
Chebyshev discrepancy

\[
 \Delta_H(q)=
 \sum_{0<|r-q|<H}
 \operatorname{sgn}(r-q)
 \left(1-\frac{|r-q|}{H}\right)\Lambda(r),
\]

where only prime powers contribute because \(\Lambda(r)=0\) otherwise. Then,
uniformly at prime-power centers \(q\to\infty\),

\[
 \boxed{
 E_h(q)=\frac{\Delta_H(q)}{2\sqrt q}
       +O_\lambda\!\left(\frac{(\log q)^3}{\sqrt q}\right).
 }
 \tag{1.1}
\]

Thus the empirical scale

\[
 E_h(q)\asymp \frac{\log q}{q^{1/4}}
\]

is exactly the scale predicted by square-root cancellation in
\(\Delta_H(q)\):

\[
 \boxed{
 |\Delta_H(q)|\ll q^{1/4}(\log q)^A.
 }
 \tag{1.2}
\]

Equation (1.1) is an unconditional reduction. Estimate (1.2), uniformly for
**every** prime-power center, is the new analytic statement that is not
currently proved.

## 2. Starting identities from the paper

The paper uses

\[
 \mathfrak P(t)=C-\mathfrak B(t),
 \qquad
 E_h=\bar p-
 \frac{\mathfrak P(t+h)-\mathfrak P(t-h)}{2h}.
\]

Its regularized residual is

\[
 \widehat{\mathfrak B}(t)
 =\mathfrak B(t)-R_0(t),
 \qquad
 R_0(t)=\sum_{k\ge1}
 \frac{e^{-(2k+1/2)t}}{2k(2k+1)},
\]

with characteristic coordinates

\[
 \widehat{\mathfrak B}(t)
 =e^{t/2}U(t)-e^{-t/2}V(t),
\]

\[
 U(t)=S_\Lambda(e^t)-t+\gamma,
 \qquad
 V(t)=\psi(e^t)-e^t+\log(2\pi).
\]

Between prime powers,

\[
 U'=-1,\qquad V'=-e^t,
\]

and at \(r=p^m\),

\[
 \Delta U=\frac{\Lambda(r)}r,
 \qquad
 \Delta V=\Lambda(r).
\]

These are equations
`Bhat-U`--`Bhat-reconstruction` in the current TeX paper.

## 3. Exact hyperbolic prime-side formula

At the central event \(q\), use half weight at the central atom:

\[
 S_\Lambda^\circ(q)
 =\sum_{r<q}\frac{\Lambda(r)}r+\frac{\Lambda(q)}{2q},
 \qquad
 \psi^\circ(q)
 =\sum_{r<q}\Lambda(r)+\frac{\Lambda(q)}2,
\]

\[
 U^\circ=S_\Lambda^\circ(q)-t+\gamma,
 \qquad
 V^\circ=\psi^\circ(q)-q+\log(2\pi),
\]

and

\[
 \widehat b^\circ(q)
 =\frac12\left(\sqrt q\,U^\circ+\frac{V^\circ}{\sqrt q}\right)
 =\frac{\widehat{\mathfrak B}'(t^-)
        +\widehat{\mathfrak B}'(t^+)}2.
\]

Define the signed hyperbolic atomic kernel

\[
 W_h(u)=
 \frac2h\,e^{-u/2}
 \sinh\!\left(\frac{h-|u|}{2}\right)
 \operatorname{sgn}(u)\,
 \mathbf 1_{\{0<|u|<h\}}.
\tag{3.1}
\]

The exact multiplicative-window discrepancy is

\[
 \Delta_{\sinh}(q,h)
 =\sum_{r=p^m}\Lambda(r)\,
 W_h\!\left(\log\frac rq\right),
\]

or, written as its two neighboring halves,

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
\tag{3.2}
\]

Also put

\[
 \kappa(h)=\frac{2\sinh(h/2)}h-1,
 \qquad
 \beta(h)=\frac{2\sinh(h/2)}h-\cosh(h/2),
\]

and, with \(\alpha_k=2k+1/2\),

\[
 \mathcal R_0(q,h)=
 \sum_{k\ge1}\frac{q^{-\alpha_k}}{2k(2k+1)}
 \left(\alpha_k-\frac{\sinh(\alpha_kh)}h\right).
\]

Then the following identity is exact:

\[
 \boxed{
 E_h(q)=
 \frac{\Delta_{\sinh}(q,h)}{2\sqrt q}
 +\kappa(h)\widehat b^\circ(q)
 +\sqrt q\,\beta(h)
 +\mathcal R_0(q,h).
 }
 \tag{3.3}
\]

### Derivation

For \(D_hf=[f(t+h)-f(t-h)]/(2h)\), the identity
\(\mathfrak P=C-\mathfrak B\) gives

\[
 E_h=D_h\mathfrak B-\mathfrak B'{}^\circ.
\]

Write the endpoint coordinates relative to the half-weight central state.
For \(r=q e^u\) on the right, its contribution to
\(\widehat{\mathfrak B}(t+h)\) is

\[
 2\frac{\Lambda(r)}{\sqrt r}
 \sinh\!\left(\frac{h-u}{2}\right).
\]

For \(r=q e^{-u}\) on the left, its contribution to
\(\widehat{\mathfrak B}(t-h)\) is the same expression and therefore enters
the centered difference with the opposite sign. The central atom contributes
equally to both endpoints and cancels, which independently verifies that its
odd triangular weight is zero.

The central state contributes

\[
 \left(\frac{2\sinh(h/2)}h-1\right)\widehat b^\circ(q),
\]

while the deterministic flows \(U'=-1\) and \(V'=-e^t\) contribute

\[
 \sqrt q\left(\frac{2\sinh(h/2)}h-\cosh(h/2)\right).
\]

Finally \(D_hR_0-R_0'\) is exactly \(\mathcal R_0(q,h)\). Combining the four
pieces gives (3.3).

This also resolves the smooth/impulse split in the paper's original
\(\mathcal T^\triangle\)-formula. There the atoms initially have linear
weights \(J_r(h-|u|)/(2h)\), while the odd integral contains
\(\mathfrak P/4+k\). Solving the regularized characteristic coordinates
recombines the \(\mathfrak P/4\) part with those linear weights, changing them
exactly into

\[
 \frac{J_r}{h}\sinh\!\left(\frac{h-|u|}{2}\right).
\]

The remaining smooth forcing is precisely the \(\sqrt q\,\beta(h)\) term;
the inherited central state and the trivial-zero correction are the other
two explicit terms in (3.3). No smooth integral remains hidden in
\(\Delta_{\sinh}\).

## 4. All non-impulse terms are already harmless

As \(h\to0\),

\[
 \kappa(h)=\frac{h^2}{24}+O(h^4),
 \qquad
 \beta(h)=-\frac{h^2}{12}+O(h^4).
\]

The classical elementary estimates

\[
 \psi(x)\ll x,
 \qquad
 \sum_{n\le x}\frac{\Lambda(n)}n=\log x+O(1)
\]

give

\[
 \widehat b^\circ(q)\ll\sqrt q.
\]

Consequently, for fixed \(\lambda\),

\[
 \kappa(h)\widehat b^\circ(q)
 +\sqrt q\,\beta(h)
 +\mathcal R_0(q,h)
 =O_\lambda\!\left(\frac{(\log q)^2}{\sqrt q}\right).
 \tag{4.1}
\]

Thus no zero estimate, RH assumption, or delicate PNT estimate is needed for
the central-state and smooth terms. The only term capable of producing the
visible cloud in the graph is the local signed prime-power discrepancy
\(\Delta_{\sinh}\).

## 5. Reduction to the ordinary additive triangular weight

Inside the window, \(|r-q|/q=O(h)\). Expanding
\(\log(r/q)\), \(\sqrt{q/r}\), and the hyperbolic sine in (3.2) shows that

\[
 \Delta_{\sinh}(q,h)=\Delta_H(q)+O_\lambda((\log q)^3).
 \tag{5.1}
\]

This estimate is elementary. The pointwise weight error is \(O_\lambda(h)\),
there are \(O(H)\) possible integer locations, and
\(\Lambda(r)\le\log(2q)\). The difference between the additive endpoints
\(q\pm H\) and the multiplicative endpoints \(qe^{\pm h}\) has length
\(O(qh^2)=O_\lambda((\log q)^2)\), giving the same total error.

Equations (3.3), (4.1), and (5.1) prove (1.1).

## 6. The exact missing deterministic estimate

On a fixed positive \(\mathfrak P\)-scale,

\[
 Q(t)^2
 =\mathfrak P(t)
 \left(2\sqrt q-\frac C2+\frac{\mathfrak P(t)}4+o(1)\right),
\]

so \(Q(t)\asymp q^{1/4}\). In that regime, (1.1) shows that

\[
 \frac{|E_h(q)|}{Q(t)}
 \ll\frac{(\log q)^A}{\sqrt q}
\]

would follow from the uniform square-root-cancellation estimate (1.2).

Without a lower bound for \(\mathfrak P(t)\), the robust state-adaptive
version must be stated directly for the exact discrepancy:

\[
 \boxed{
 |E_h(q)|\ll Q(t)\frac{(\log q)^A}{\sqrt q}
 }
 \tag{6.1}
\]

The additive statement
\(|\Delta_H(q)|\ll Q(t)(\log q)^A\) implies (6.1) only with enough room to
absorb the explicit error in (1.1), or with a separate relative bound for
that error. This qualification is essential. Under RH, recurrence gives
small values of \(\mathfrak P(t)\) on the continuous axis, but it does not
place them at \(t=\log(p^a)\). In either case, no fixed lower bound
\(Q\gg q^{1/4}\) is available for a universal eventwise proof.

For preservation of the lower cone, an absolute-value estimate is stronger
than necessary. It would be enough to prove the corresponding one-sided
lower-tail estimate for \(\Delta_H(q)\).

The logical status is therefore:

- (3.3), (4.1), (5.1), and hence (1.1), are unconditional;
- converting an absolute \(E_h\)-bound into the plotted \(E_h/Q\)-bound
  requires a quantitative lower bound for \(Q(t)\);
- replacing \(Q(t)\) by \(q^{1/4}\) assumes that
  \(\mathfrak P(t)\) is bounded below on the events under consideration;
- neither (1.2) nor the exact state-adaptive target (6.1) has been proved.

## 7. Why simple absolute bounds do not prove the observed law

The two unsigned halves each have weighted von Mangoldt mass of order

\[
 H=\lambda\sqrt q\log q.
\]

Brun--Titchmarsh for the prime contribution, plus an elementary count of
higher prime powers, gives only

\[
 |\Delta_H(q)|\ll H,
 \qquad
 |E_h(q)|\ll\log q.
 \tag{7.1}
\]

This loses a factor of approximately \(q^{1/4}\) relative to (1.2). Bounding
the right and left halves separately therefore cannot prove the slope seen in
the graph. It can give
\(|E_h|/Q\ll (\log q)/q^{1/4}\) only on a subsequence where an independent
lower bound \(Q\gg q^{1/4}\) is already known; such a lower bound is not
available uniformly.

The desired factor comes from cancellation **between** the two adjacent
weighted halves, not from a sharper upper bound for either half.

## 8. Why the numerical slope is nevertheless natural

The window contains about

\[
 \frac{H}{\log q}\asymp\sqrt q
\]

primes. A prime contributes weight of size \(\log q\). Square-root
cancellation therefore predicts

\[
 |\Delta_H(q)|_{\mathrm{typ}}
 \asymp \sqrt{\sqrt q}\,\log q
 =q^{1/4}\log q.
\]

Substitution into (1.1) gives

\[
 |E_h(q)|_{\mathrm{typ}}\asymp\frac{\log q}{q^{1/4}},
 \qquad
 \frac{|E_h(q)|}{Q(t)}\asymp\frac{\log q}{\sqrt q}
\]

whenever \(Q(t)\asymp q^{1/4}\). This explains the downward-moving point
cloud without treating a finite regression as a proof.

This scale is consistent with the short-interval variance picture of
Montgomery--Soundararajan, but their distributional/average framework does
not supply a uniform bound at every prime-power center:

- Hugh Montgomery and K. Soundararajan,
  [Primes in short intervals](https://arxiv.org/abs/math/0409258).
- Kevin Ford, K. Soundararajan, and Alexandru Zaharescu,
  [On the distribution of imaginary parts of zeros of the Riemann zeta
  function, II](https://arxiv.org/abs/0805.2745).

## 9. Numerical identity check

The reproducible regression test is
[test_prime_side_eh_decomposition.py](test_prime_side_eh_decomposition.py).
Run it from this companion directory with

    python -m unittest test_prime_side_eh_decomposition.py -v

It checks (3.3) through two independent binary64 evaluation paths at the
prime and higher-prime-power centers
\(q\in\{2,5,16,125,997,4096,7919\}\) and
\(\lambda\in\{1/4,1/2,1\}\).

As an independent binary64 check, at \(q=997\) and \(\lambda=1/2\), the direct
prime-side endpoint formula gives

\[
 E_h=0.059361592425059694.
\]

The four terms in (3.3) are respectively

\[
 0.09086018560044935,\quad
 -0.00003292259346748,\quad
 -0.03146567041605164,\quad
 -1.6593713401\times10^{-10},
\]

whose sum is \(0.05936159242499311\). The difference from the independent
endpoint calculation is \(6.7\times10^{-14}\), consistent with binary64
roundoff.

## 10. Conclusion

The finite graph and full-population audit identify a coherent scaling
candidate. The new exact reduction is

\[
 E_h(q)=\frac{1}{2\sqrt q}
 \left[\text{right weighted prime mass}
       -\text{left weighted prime mass}\right]
 +O_\lambda\!\left(\frac{(\log q)^3}{\sqrt q}\right).
\]

What remains in this discrepancy branch is not the smooth integral: that
part is now bounded elementarily. The arithmetic obstacle is a uniform,
signed, square-root-cancellation estimate for the two neighboring triangular
prime masses, with enough relative control to imply the exact
state-adaptive target (6.1). Current absolute short-interval estimates do not
contain this cancellation. The coupled central-secant branch remains open as
well.
