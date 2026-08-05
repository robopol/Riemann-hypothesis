# Prime-power cell dynamics of the residual function \(B(t)\)

**Date:** 2026-08-04  
**Status:** research note, not a proof of the Riemann hypothesis  
**Scope:** exact hybrid dynamics, polynomial-energy search, and exploratory
binary64 tests. The two original RH/MVDC papers were not edited.

## 1. Main outcome

The function

\[
B(t)=\sum_\rho
\frac{e^{(\rho-1/2)t}}{\rho(1-\rho)}
\]

can be analyzed directly from the prime side without generating zeros or CA
profiles. The resulting system is especially simple.

1. Between consecutive prime powers, \(B\) satisfies an explicit linear ODE.
   At a prime power \(q=p^m\), \(B\) remains continuous and its derivative
   receives the exact positive impulse

   \[
   \Delta B'(\log q)=\frac{\Lambda(q)}{\sqrt q}
   =\frac{\log p}{\sqrt q}.
   \]

2. After removing an explicit exponentially decaying trivial-zero term, the
   dynamics becomes

   \[
   \widehat B''-\frac14\widehat B
   =-e^{t/2}
   +\sum_q\frac{\Lambda(q)}{\sqrt q}\,
   \delta(t-\log q).
   \]

   This is an exact forced hyperbolic oscillator whose event schedule consists
   of the prime powers.

3. The natural two-dimensional state is exactly the pair of classical prime
   errors

   \[
   U(t)=\sum_{q\le e^t}\frac{\Lambda(q)}q-t+\gamma_E,
   \qquad
   V(t)=\psi(e^t)-e^t+\log(2\pi).
   \]

   Thus the cell dynamics has not removed the prime fluctuation. It has located
   the precise cancellation between the two prime errors that a successful
   proof must control.

4. No time-independent coercive quadratic energy can be decreasing on every
   cell and at every impulse. A polynomially degenerating family does make the
   unavoidable square-jump budget summable, but a signed prime-slope work term
   remains.

5. A more promising one-sided formulation is obtained from

   \[
   \Phi(t)=C-B(t),
   \qquad
   C=B(0)=2+\gamma_E-\log(4\pi).
   \]

   In fact, even eventual non-negativity of \(\Phi\) is equivalent to RH. For
   \(t\ge\log 2\), define

   \[
   R(t)=\sqrt{\frac{\Phi(t)^2}{4}+2k(t)\Phi(t)},
   \qquad
   H(t)=\Phi'(t)^2-R(t)^2.
   \]

   While \(\Phi>0\), \(H\) strictly decreases inside every open cell. At a
   prime power, preservation of the cone \(H<0\) reduces to the single exact
   kick inequality

   \[
   \boxed{
   R(\log q)+\Phi'(\log q^-)
   >\frac{\Lambda(q)}{\sqrt q}.
   }
   \]

   Proving this inequality for every prime power, together with the explicit
   initial state at \(q=2\), would prove \(\Phi(t)>0\) eventually and hence RH.

6. An exploratory event-driven scan processed all \(3\,342\,115\) prime-power
   events through \(56\,048\,351\), as well as \(23\,746\) stationary points
   inside the cells. The kick inequality and the cone \(H<0\) held throughout
   the scan. This is finite binary64 evidence, not an interval certificate and
   not a proof of the universal inequality.

The concrete advance is therefore a reduction, not a completed RH proof: the
prime-power dynamics and a falsifiable discrete invariant are now explicit.
The remaining task is to prove the kick inequality, or a block-averaged version
that implies it cannot fail.

## 2. Starting formula

Let

\[
S_\Lambda(x)=\sum_{n\le x}\frac{\Lambda(n)}n,
\qquad
\psi(x)=\sum_{n\le x}\Lambda(n),
\]

and put

\[
T(z)=z\log(1-z^2)-2z+\log\frac{1+z}{1-z}.
\]

The exact prime-side formula established in the signed-cancellation study is

\[
\boxed{
B(\log x)=\sqrt{x}\left[
S_\Lambda(x)-\log x+\gamma_E
-\left(\frac{\psi(x)}x-1\right)
-\frac{\log(2\pi)}x
-\frac12T(1/x)
\right].
}
\]

The two jumps at a prime power cancel in the value of \(B\):

\[
\Delta S_\Lambda(q)=\frac{\Lambda(q)}q,
\qquad
\Delta\frac{\psi(q)}q=\frac{\Lambda(q)}q.
\]

This continuity is what makes an exact cell analysis possible.

## 3. Exact open-cell formula

Write

\[
t=\log x,
\qquad z=e^{-t},
\qquad c=\log(2\pi).
\]

On an open cell between consecutive prime powers, the quantities

\[
S=S_\Lambda(e^t),
\qquad P=\psi(e^t)
\]

are constant. Set

\[
a=\gamma_E+1+S,
\qquad d=c+P.
\]

Then

\[
\boxed{
B(t)=e^{t/2}(a-t)-d e^{-t/2}
-\frac12e^{t/2}T(e^{-t}).
}
\]

The elementary identities

\[
T'(z)=\log(1-z^2),
\qquad
T''(z)=-\frac{2z}{1-z^2}
\]

give

\[
\boxed{
\begin{aligned}
B'(t)={}&e^{t/2}\left[
\frac{a-t}{2}-1-\frac{T(z)}4
+\frac z2\log(1-z^2)
\right]
+\frac d2e^{-t/2}.
\end{aligned}
}
\]

A second differentiation yields the exact cell equation

\[
\boxed{
B''(t)-\frac14B(t)
=-g(t),
}
\]

where

\[
g(t)=e^{t/2}\left(
1-\frac{e^{-3t}}{1-e^{-2t}}
\right)
=e^{t/2}-\frac{e^{-5t/2}}{1-e^{-2t}}.
\]

For \(t\ge\log2\), both \(g(t)>0\) and \(g'(t)>0\).

## 4. Exact prime-power event law

Let \(q=p^m\), \(\tau_q=\log q\), and

\[
J_q=\frac{\Lambda(q)}{\sqrt q}
=\frac{\log p}{\sqrt q}.
\]

At \(q\), the cell constants jump by

\[
\Delta a=\frac{\Lambda(q)}q,
\qquad
\Delta d=\Lambda(q).
\]

Substitution in the formulas above gives

\[
\boxed{
\Delta B(\tau_q)=0,
\qquad
\Delta B'(\tau_q)=J_q>0.
}
\]

Consequently, in the distributional sense on \((0,\infty)\),

\[
\boxed{
B''-\frac14B
=-e^{t/2}+\frac{e^{-5t/2}}{1-e^{-2t}}
+\sum_{q=p^m}J_q\,\delta(t-\tau_q).
}
\]

The ordinary left and right values of \(B''\) agree at an event. The singular
part is exactly the delta mass caused by the jump of \(B'\).

### 4.1 Exact impulse representation

Define the no-prime archimedean trajectory

\[
A_0(t)=e^{t/2}(\gamma_E+1-t)
-c e^{-t/2}
-\frac12e^{t/2}T(e^{-t}).
\]

The complete solution is

\[
\boxed{
B(t)=A_0(t)
+2\sum_{\tau_q\le t}J_q
\sinh\frac{t-\tau_q}{2}.
}
\]

With the right derivative convention,

\[
\boxed{
B'(t+)=A_0'(t)
+\sum_{\tau_q\le t}J_q
\cosh\frac{t-\tau_q}{2}.
}
\]

These formulas reproduce both continuity and the derivative impulse exactly.

### 4.2 The origin is singular in the derivative

Although

\[
B(0)=C=2+\gamma_E-\log(4\pi),
\]

one must not initialize the state as \((B(0),0)\). In fact,

\[
\boxed{
B(t)=C+\frac t2\log t
+\frac{\gamma_E+\log(2\pi)-1}{2}\,t+o(t),
\qquad t\downarrow0,
}
\]

and therefore \(B'(0+)=-\infty\). The event dynamics below is started at
\(t=\log2\), where both one-sided derivatives are finite.

## 5. Removing the trivial-zero forcing

Define

\[
\boxed{
\widehat B(t)=B(t)+\frac12e^{t/2}T(e^{-t}).
}
\]

Then

\[
\boxed{
\widehat B(t)
=e^{t/2}\bigl(\gamma_E+1+S_\Lambda(e^t)-t\bigr)
-e^{-t/2}\bigl(c+\psi(e^t)\bigr),
}
\]

and its hybrid equation is particularly clean:

\[
\boxed{
\widehat B''-\frac14\widehat B
=-e^{t/2}+\sum_qJ_q\delta(t-\tau_q).
}
\]

Unlike \(B'\), the regularized derivative has finite initial data:

\[
\widehat B(0)=\gamma_E+1-c,
\qquad
\widehat B'(0)=\frac{\gamma_E+c-1}{2}.
\]

The removed term is explicit. Since \(T(0)=0\) and
\(T'(z)=\log(1-z^2)\),

\[
T(z)=-\sum_{k\ge1}\frac{z^{2k+1}}{k(2k+1)}.
\]

Thus

\[
\boxed{
B(t)-\widehat B(t)
=\sum_{k\ge1}
\frac{e^{-(2k+1/2)t}}{2k(2k+1)}
=\frac16e^{-5t/2}+O(e^{-9t/2}).
}
\]

It is positive and exponentially decaying. Hence \(B\) and \(\widehat B\)
have exactly the same polynomial or subexponential growth question.

## 6. Exact two-dimensional prime-error state

Put

\[
\boxed{
U(t)=e^{-t/2}\left(\widehat B'(t)+\frac12\widehat B(t)\right)
=S_\Lambda(e^t)-t+\gamma_E,
}
\]

and

\[
\boxed{
V(t)=e^{t/2}\left(\widehat B'(t)-\frac12\widehat B(t)\right)
=\psi(e^t)-e^t+c.
}
\]

The reconstruction formulas are

\[
\boxed{
\widehat B=e^{t/2}U-e^{-t/2}V,
\qquad
\widehat B'=\frac12\left(e^{t/2}U+e^{-t/2}V\right).
}
\]

On every open cell,

\[
U'=-1,
\qquad
V'=-e^t,
\]

and at \(q\),

\[
\Delta U=\frac{\Lambda(q)}q,
\qquad
\Delta V=\Lambda(q).
\]

For consecutive prime powers \(q_k<q_{k+1}\), the exact right-state event
map is therefore

\[
\boxed{
U_{k+1}=U_k-\log\frac{q_{k+1}}{q_k}
+\frac{\Lambda(q_{k+1})}{q_{k+1}},
}
\]

\[
\boxed{
V_{k+1}=V_k-(q_{k+1}-q_k)
+\Lambda(q_{k+1}).
}
\]

This is the requested exact prime-power dynamics. It contains no zero
computation. It also explains why the problem remains hard: polynomial control
of \(\widehat B\) requires exponentially accurate cancellation between the two
normalized prime errors.

### 6.1 Full-\(B\) coordinates

For completeness, without removing the trivial term one obtains

\[
e^{-t/2}\left(B'+\frac B2\right)
=S_\Lambda(e^t)-t+\gamma_E+e^{-t}-\operatorname{atanh}(e^{-t}),
\]

\[
e^{t/2}\left(B'-\frac B2\right)
=\psi(e^t)-e^t+c+\frac12\log(1-e^{-2t}).
\]

The same event jumps occur; only the explicit smooth forcing changes.

## 7. A nonlocal identity exposing the cancellation

Let

\[
E(x)=\psi(x)-x.
\]

Partial summation followed by the PNT boundary condition gives

\[
S_\Lambda(x)-\log x+\gamma_E
=\frac{E(x)}x-\int_x^\infty\frac{E(u)}{u^2}\,du.
\]

The endpoint \(E(x)/x\) cancels in \(\widehat B\), leaving

\[
\boxed{
\widehat B(\log x)
=-\sqrt{x}\int_x^\infty
\frac{\psi(u)-u}{u^2}\,du
-\frac{\log(2\pi)}{\sqrt{x}}.
}
\]

This formula is conceptually important. The smallness of \(B\) is a future
signed cancellation, not a pointwise smallness statement for \(\psi(x)-x\).
A forward local energy must somehow encode this terminal cancellation.

In particular, the target

\[
\widehat B(\log x)=O((\log x)^A)
\]

is equivalent to the RH-strength tail estimate

\[
\int_x^\infty\frac{\psi(u)-u}{u^2}\,du
=O\left(x^{-1/2}(\log x)^A\right).
\]

Taking absolute values inside the integral destroys the cancellation and does
not reach this scale.

## 8. Why a polynomial envelope would settle RH

The preceding signed-cancellation note proved

\[
\mathrm{RH}
\iff B\text{ is bounded on }\mathbb R.
\]

The same Laplace-transform argument gives the stronger practical statement:
any subexponential estimate

\[
|B(t)|\le e^{o(t)}
\]

already excludes a zero with \(\Re\rho>1/2\). In particular,

\[
\boxed{
|B(t)|\le K(1+t)^A
\quad(t\ge0)
\quad\Longrightarrow\quad\mathrm{RH}.
}
\]

Because \(B-\widehat B\) decays exponentially, the same is true with
\(\widehat B\) in place of \(B\).

Thus the proposed polynomial-energy route is logically valid. Its difficulty
is not choosing a sufficiently large exponent \(A\); it is proving any global
subexponential bound from the prime-power dynamics.

## 9. A stronger one-sided equivalent: \(\Phi=C-B\)

Define

\[
\Phi(t)=C-B(t).
\]

Under RH, pairing ordinates \(\gamma\) and \(-\gamma\) gives

\[
\Phi(t)
=\sum_{\gamma>0}
\frac{2(1-\cos\gamma t)}{\gamma^2+1/4}
\ge0.
\]

The converse also holds, and eventual positivity is enough:

\[
\boxed{
\mathrm{RH}
\iff \Phi(t)\ge0\text{ for all sufficiently large }t.
}
\]

Here is the Laplace-transform proof of the converse. For \(\Re s>1/2\),

\[
F(s)=\int_0^\infty e^{-st}\Phi(t)\,dt
=\frac Cs
-\sum_\rho
\frac{1}{\rho(1-\rho)[s-(\rho-1/2)]}.
\]

The series on the right is meromorphic and normally convergent off its poles;
on compact sets its high summands are \(O(|\gamma|^{-3})\). Suppose that
\(\Phi(t)\ge0\) for \(t\ge T\), and apply the classical Landau theorem for
Laplace transforms of nonnegative functions to

\[
G(s)=\int_0^\infty e^{-su}\Phi(T+u)\,du.
\]

If its abscissa of convergence were a positive number \(\sigma_c\), the real
point \(s=\sigma_c\) would be singular. The displayed meromorphic continuation
is analytic at every positive real \(s\), because \(\zeta\) has no real zero in
\((0,1)\). Hence \(\sigma_c\le0\), so \(G\) is holomorphic in \(\Re s>0\).
The identity theorem now excludes every pole \(s=\rho-1/2\) with positive real
part. The functional equation excludes the reflected zeros with real part
below \(1/2\).

This one-sided criterion belongs to the same general Weil/screw-function
framework as Suzuki's criterion. Suzuki's function uses the weight
\(1/\gamma^2\), whereas the present \(\Phi\) uses
\(1/(\gamma^2+1/4)\). The functions are related but not identical, and no
novelty claim is made here for the positivity principle.

## 10. Cell geometry and exact extrema

### 10.1 Geometry of \(\widehat B\)

On a fixed cell, \(\widehat B'=0\) is equivalent to

\[
e^t(t+2-a)=d.
\]

Therefore the cell contains at most one stationary point, explicitly

\[
\boxed{
t_*=a-2+W_0\left(de^{2-a}\right),
}
\]

where \(W_0\) is the principal Lambert function. Let

\[
w=W_0(de^{2-a}).
\]

If \(t_*\) lies in the cell, then

\[
\widehat B(t_*)=2e^{t_*/2}(1-w),
\]

and

\[
\boxed{
\widehat B''(t_*)
=-\frac12e^{t_*/2}(1+w)<0.
}
\]

Thus every interior stationary point of \(\widehat B\) is a strict maximum.
There are no interior minima. Since every event jumps the derivative upward,
an event can create a local minimum but never a local maximum.

This reduces an exact envelope scan to:

- the continuous values at prime-power events;
- at most one Lambert-\(W\) maximum inside each cell.

The full function \(B\) differs only by the explicit positive trivial-zero
correction. The exploratory scan below locates stationary points of the full
\(B\), rather than substituting this asymptotic equivalence.

### 10.2 Uniqueness of the full-\(B\) cell maximum

The full function also has at most one stationary point in every cell with
\(t\ge\log2\). Define

\[
h(z)=z+\frac z2\log(1-z^2)
-\frac12\log\frac{1+z}{1-z}
\]

and

\[
F(t)=t+2-a-de^{-t}-h(e^{-t}).
\]

The derivative formula is exactly

\[
B'(t)=-\frac12e^{t/2}F(t).
\]

Moreover,

\[
\frac{d}{dt}h(e^{-t})
=-\frac z2\log(1-z^2)
+\frac{2z^3}{1-z^2}.
\]

For \(0<z\le1/2\), the elementary inequality

\[
-\log(1-z^2)\le\frac{z^2}{1-z^2}
\]

gives

\[
0<\frac{d}{dt}h(e^{-t})
\le\frac{5z^3}{2(1-z^2)}
\le\frac5{12}.
\]

Since \(d>0\),

\[
\boxed{
F'(t)=1+de^{-t}-\frac{d}{dt}h(e^{-t})>0.
}
\]

Thus \(B'\) can vanish at most once in the cell, and any zero is a strict
positive-to-negative crossing because

\[
B''(t_*)=-\frac12e^{t_*/2}F'(t_*)<0.
\]

Consequently, from \(q=2\) onward the full \(B\), not only \(\widehat B\), has
no interior cell minimum. An endpoint-sign change identifies the unique cell
maximum, so the event-driven root scan is exhaustive on the scanned cells.

## 11. Quadratic energy audit

### 11.1 Why a fixed coercive energy cannot work locally

The homogeneous cell flow for \((\widehat B,\widehat B')\) has modes
\(e^{t/2}\) and \(e^{-t/2}\). Over a cell of length \(h\), its eigenvalues are
\(e^{h/2}\) and \(e^{-h/2}\). Along the unstable eigenvector, every constant
positive-definite quadratic form grows by \(e^h\). Hence no such form can be
monotone decreasing on every possible open-cell state.

At an event, the map is

\[
(\widehat B,\widehat B')
\longmapsto
(\widehat B,\widehat B'+J_q).
\]

Any quadratic form with a positive \(\widehat B'^2\) coefficient has a jump
containing a positive \(J_q^2\) term and an indefinite state-linear term. A
successful energy must therefore use the actual correlation between the prime
gaps and the state, become time-dependent or nearly degenerate, or be nonlocal.

### 11.2 The normalized-mode family

Consider

\[
E_\eta=e^tU^2+e^{-t}V^2-2\eta UV.
\]

For \(|\eta|<1\), this is coercive in the two normalized modes. At an event,

\[
\boxed{
\Delta E_\eta
=2(1-\eta)\left[
\Lambda(q)\left(U+\frac Vq\right)
+\frac{\Lambda(q)^2}{q}
\right].
}
\]

The unique jump-invariant member is \(\eta=1\), for which

\[
E_1=\widehat B^2.
\]

It loses derivative coercivity. This degeneracy is structural: in the variables

\[
\alpha=e^{t/2}U,
\qquad
\beta=e^{-t/2}V,
\]

every event is the common translation

\[
(\alpha,\beta)\mapsto(\alpha+J_q,\beta+J_q).
\]

A positive semidefinite quadratic form invariant under every such translation
must annihilate the direction \((1,1)\), so it can control only
\(\alpha-\beta=\widehat B\), not the derivative.

### 11.3 A polynomially degenerating energy

Put

\[
\eta(t)=1-\varepsilon(t),
\qquad
\varepsilon(t)=(1+t)^{-r}.
\]

Then

\[
\boxed{
E_\varepsilon
=\left(1-\frac\varepsilon2\right)\widehat B^2
+2\varepsilon\widehat B'^2.
}
\]

It remains coercive in \(\widehat B\), while derivative coercivity is allowed
to weaken only polynomially. On an open cell,

\[
\boxed{
E_\varepsilon'
=2\widehat B\widehat B'
+\varepsilon'\left(
2\widehat B'^2-\frac12\widehat B^2
\right)
-4\varepsilon e^{t/2}\widehat B'.
}
\]

At an event,

\[
\boxed{
\Delta E_\varepsilon
=4\varepsilon J_q\widehat B'_-
+2\varepsilon J_q^2.
}
\]

The purely positive square-jump budget is

\[
\sum_{q=p^m}
(1+\log q)^{-r}\frac{\Lambda(q)^2}{q}.
\]

The terms with \(m\ge2\) converge for every fixed \(r\). By the PNT and
partial summation, the prime layer has the convergence behavior

\[
\int^\infty v^{1-r}\,dv.
\]

Consequently,

\[
\boxed{
\sum_q(1+\log q)^{-r}
\frac{\Lambda(q)^2}{q}<\infty
\iff r>2.
}
\]

This is a real gain: the unavoidable positive quadratic kick cost can be made
finite. The unresolved term is the signed work

\[
4\varepsilon J_q\widehat B'_-.
\]

Controlling its cumulative interaction with the open-cell drift is the exact
remaining problem for this polynomial-energy family.

## 12. The first-crossing cone

The one-sided variable \(\Phi=C-B\) satisfies, on an open cell,

\[
\boxed{
\Phi''-\frac14\Phi=k(t),
}
\]

where

\[
k(t)=g(t)-\frac C4
=e^{t/2}
-\frac{e^{-t/2}}{e^{2t}-1}
-\frac C4.
\]

For \(t\ge\log2\),

\[
k(t)>0,
\qquad
k'(t)>0.
\]

At \(q=p^m\),

\[
\Delta\Phi=0,
\qquad
\Delta\Phi'=-J_q.
\]

Define

\[
\boxed{
H=\Phi'^2-\frac14\Phi^2-2k\Phi.
}
\]

Direct differentiation gives the exact smooth-cell law

\[
\boxed{
H'=-2k'(t)\Phi(t).
}
\]

At an event,

\[
\boxed{
\Delta H=J_q^2-2J_q\Phi'(\tau_q^-).
}
\]

The event increments do not have a fixed sign. The useful invariant is not
monotonicity at every event, but the cone \(H<0\).

Assume \(\Phi>0\), and put

\[
R=\sqrt{\frac{\Phi^2}{4}+2k\Phi}.
\]

Then

\[
H<0\iff |\Phi'|<R.
\]

If \(H<0\) immediately after an event, a first zero of \(\Phi\) inside the
next cell is impossible: until that putative first zero, \(H\) decreases, but
at the zero one would have \(H=\Phi'^2\ge0\).

At the next event, let \(v_-=\Phi'(\tau_q^-)\). The pre-event cone gives
\(-R<v_-<R\), while the kick gives \(v_+=v_--J_q\). The upper cone boundary is
automatic. The lower boundary is preserved exactly when

\[
\boxed{
R+v_->J_q.
}
\]

Therefore the following is a complete discrete sufficient condition:

1. verify \(\Phi(\log2)>0\) and \(H(\log2+)<0\);
2. prove

   \[
   R(\log q)+\Phi'(\log q^-)
   >\frac{\Lambda(q)}{\sqrt q}
   \]

   at every later prime power.

Induction then gives \(\Phi(t)>0\) for every \(t\ge\log2\). By the eventual
positivity criterion, this proves RH. The singular derivative at \(t=0\) is
irrelevant to this induction.

### 12.1 Scalar reserve dynamics

The cone can be expressed as a particularly transparent reserve-withdrawal
system. Put

\[
m=R+\Phi',
\qquad
n_*=R-\Phi',
\qquad
-H=mn_*.
\]

On an open cell, let

\[
A=k+\frac\Phi4.
\]

Differentiating \(R^2=\Phi^2/4+2k\Phi\) gives the exact equations

\[
\boxed{
m'=\frac ARm+\frac{k'\Phi}{R}>0,
}
\]

and

\[
\boxed{
n_*'=-\frac ARn_*+\frac{k'\Phi}{R}.
}
\]

At a prime-power event,

\[
\boxed{
m^+=m^--J_q,
\qquad
n_*^+=n_*^-+J_q.
}
\]

Thus the smooth flow strictly replenishes the lower-bound reserve \(m\), while
the event withdraws exactly \(\Lambda(q)/\sqrt q\). The universal target is
the scalar inequality

\[
\boxed{
m(\log q^-)>J_q.
}
\]

This is equivalent to the kick condition above and is the cleanest state
variable for future block experiments.

At \(q=2\), the binary64 values are

\[
\begin{aligned}
\Phi&\approx0.06352355985,\\
\Phi'_-&\approx0.22501526882,\\
J_2&\approx0.49012907173,\\
R&\approx0.38635237215,\\
\Phi'_+&\approx-0.26511380292,\\
R+\Phi'_+&\approx0.12123856924,\\
H_+&\approx-0.07898282697.
\end{aligned}
\]

These initial values can be certified separately by elementary interval
arithmetic if the universal kick inequality is eventually proved.

### 12.2 Exact arithmetic form of the kick

For an event \(x=q\), put

\[
\ell(x)=\frac12\log(1-x^{-2}).
\]

The full prime-side coordinates give

\[
\boxed{
\Phi'(\log x^+)
=\frac12(\Phi-C)
+\frac{x-\psi(x)-c-\ell(x)}{\sqrt x}.
}
\]

Hence the post-event condition \(R+\Phi'^+>0\) is exactly

\[
\boxed{
\psi(x)-x
<\sqrt x\left[R+\frac{\Phi-C}{2}\right]
-c-\ell(x).
}
\]

If \(\Phi\) were already known to remain comparable to a positive constant,
the leading right-hand scale would be \(\sqrt{2\Phi}\,x^{3/4}\). Current
unconditional pointwise PNT envelopes do not provide this fixed
\(x^{3/4}\)-scale one-sided estimate. Replacing the adaptive state
\((\Phi,\Phi')\) by a crude absolute envelope is therefore a dead end; a proof
must preserve the state correlation or group the impulses into blocks.

## 13. Exploratory computation

### 13.1 Reproducibility

The event-driven implementation is

`current_support/robin_mvdc_status_support/exploratory_b_cell_dynamics.py`.

Its full binary64 report is

`current_support/robin_mvdc_status_support/exploratory_b_cell_dynamics_report.json`.

Run from the support directory with

```powershell
python exploratory_b_cell_dynamics.py --limit 56048351
```

The script:

- sieves all primes and generates every prime power through the limit;
- evaluates the exact left and right cell states;
- checks continuity and the derivative impulse law;
- compares sampled values with the independently implemented prime-side
  formula for \(B\);
- locates every derivative zero inside a cell by safeguarded bisection;
- tests polynomial normalizations, simple quadratic energies, the
  polynomially degenerating family, and the first-crossing cone.

All numerical results in this section are ordinary binary64 diagnostics. They
are not directed-rounding interval enclosures.

### 13.2 Scan size and internal checks

The full run used

\[
2\le q\le56\,048\,351,
\qquad
\log2\le t\le17.8417252869.
\]

It processed

- \(3\,341\,013\) primes;
- \(3\,342\,115\) prime-power events;
- \(23\,746\) cells containing a stationary point of the full \(B\).

The principal validation residuals were

\[
\max_q|B(q+)-B(q-)|
\approx1.46\cdot10^{-11},
\]

\[
\max_q\left|
\Delta B'-\frac{\Lambda(q)}{\sqrt q}
\right|
\approx7.50\cdot10^{-12},
\]

and the maximum discrepancy against sampled calls to the existing prime-side
implementation of \(B\) was approximately

\[
1.75\cdot10^{-11}.
\]

All report consistency checks passed.

### 13.3 Observed extrema

Over all events and all located stationary points,

\[
\min B\approx-0.0226812404648
\quad\text{at }q=319\,439,
\]

and

\[
\max B\approx0.0246942372460
\]

at the interior point

\[
x=e^t\approx34\,186\,454.4254.
\]

Thus the observed maximum absolute value was

\[
0.0246942372460<C\approx0.0461914179322.
\]

Every one of the \(23\,746\) derivative crossings was from positive to
negative; none was from negative to positive. This is consistent with the
cell-maximum geometry.

These extrema establish nothing beyond the finite scan. Their role is to test
the exact formulas and to identify the tight transitions for candidate
invariants.

### 13.4 First-crossing energy

At all scanned event states,

\[
\Phi>0
\]

and both one-sided states satisfied \(H<0\). The largest observed one-sided
value of \(H\) was

\[
H\approx-0.07898282697
\quad\text{at }q=2,
\]

while the smallest was approximately

\[
-810.70866.
\]

The endpoint value was approximately

\[
H(\log56\,048\,351+)\approx-653.50380.
\]

The smallest right-state cone margin

\[
R+\Phi'_+
\]

was approximately

\[
0.1212385692
\quad\text{at }q=2.
\]

The largest observed ratio \(|\Phi'|/R\) was approximately

\[
0.728664
\quad\text{at }q=5.
\]

Equivalently, every scanned state remained strictly inside the cone.

This happened despite the fact that \(1\,694\,123\) impulses, about
\(50.69\%\), increased \(H\). The positive event increments summed to about
\(2714.163\), the negative increments to about \(-2679.288\), and their net
was therefore positive. The open-cell drift summed to about \(-688.280\) and
had no sign violation above the numerical tolerance. The finite safety comes
from the balance of the cell drift and the impulses, not from eventwise
monotonicity.

### 13.5 Simple quadratic candidates

The scan tested fixed energies of the forms

\[
B^2+aB'^2,
\qquad
B^2+a(B'\pm B)^2,
\qquad
B^2+\frac{B'^2}{(1+t)^{2r}},
\]

over several values of \(a\) and \(r\). None was monotone across both cells
and events. Roughly half of the transitions increased each simple candidate,
as predicted by the structural no-go argument.

Finite failure rules out those particular monotonicity claims. Finite success
of another candidate would still require a universal analytic proof.

### 13.6 Polynomially degenerating energy candidates

The exact regularized family

\[
E_r(t)=\left(1-\frac{(1+t)^{-r}}2\right)\widehat B(t)^2
+2(1+t)^{-r}\widehat B'(t)^2
\]

was tested at every one-sided event state for

\[
r\in\{0,1,2,3,4,6,8\}.
\]

For every tested \(r\ge1\), the largest **event-state** energy was already at
\(q=2\). Interior maxima of this energy were not searched, so this statement
must not be read as a continuous-cell envelope. The endpoint energies were:

| \(r\) | endpoint \(E_r\) | signed jump sum | signed open-cell drift |
|---:|---:|---:|---:|
| 1 | \(6.17377\times10^{-3}\) | \(+3.92714\) | \(-3.94560\) |
| 2 | \(3.33785\times10^{-4}\) | \(+0.301676\) | \(-0.316886\) |
| 3 | \(2.38354\times10^{-5}\) | \(+0.0660037\) | \(-0.0761561\) |
| 4 | \(7.38520\times10^{-6}\) | \(+0.0313229\) | \(-0.0383213\) |
| 6 | \(6.46579\times10^{-6}\) | \(+0.00991504\) | \(-0.0139360\) |
| 8 | \(6.46320\times10^{-6}\) | \(+0.00330353\) | \(-0.00628547\) |

Thus the total finite balance was negative for every tested \(r\ge1\). For
\(r=0\), it was positive and the largest event-state energy was about
\(1.144739\) at \(q=1423\).

The favorable total balance does not supply a local Lyapunov law: roughly half
the individual event and event-to-event transitions still increased the
energy. For the proof-relevant range \(r>2\), the experiment supports a
block-cancellation search, not eventwise monotonicity.

## 14. What has and has not been achieved

### Exact advances

The following are exact, unconditional identities:

1. the open-cell ODE for \(B\) and \(\widehat B\);
2. continuity of \(B\) and the impulse
   \(\Delta B'=\Lambda(q)/\sqrt q\);
3. the two-dimensional event map in \((U,V)\);
4. the impulse representation and the future-tail identity;
5. the Lambert-\(W\) formula for the unique cell maximum of
   \(\widehat B\);
6. convergence of the square-jump budget for the polynomial energy exactly
   when \(r>2\);
7. the first-crossing energy law and its reduction to one kick inequality;
8. eventual non-negativity of \(C-B\) as an RH-equivalent target.

### What remains open

No polynomial envelope has been proved. In particular, the finite computation
does not prove that

\[
R(\log q)+\Phi'(\log q^-)
>\frac{\Lambda(q)}{\sqrt q}
\]

continues forever. Nor does it prove a global bound for the signed work in the
polynomially degenerating energy.

The exact obstruction is now sharper than “prime fluctuations are difficult”:

- the positive square-kick budget is summable for \(r>2\);
- the signed cross-work \(4\varepsilon J_q\widehat B'_-\) is not yet bounded;
- for the one-sided cone, one must prevent a downward derivative kick from
  crossing the moving lower boundary \(-R\).

Either of the last two controls would be a genuinely RH-strength new lemma.

## 15. Recommended next experiments and proof targets

### 15.1 Block the kick inequality

Instead of requiring every event increment of \(H\) to be negative, sum the
exact balance over a block \([\log q_j,\log q_k]\):

\[
H_k^+-H_j^+
=-2\int_{\log q_j}^{\log q_k}k'(t)\Phi(t)\,dt
+\sum_{j<\ell\le k}
\left(J_{q_\ell}^2
-2J_{q_\ell}\Phi'(\tau_{q_\ell}^-)
\right).
\]

The computation shows that the event sum may be locally positive while the
cell integral dominates on longer blocks. A useful next lemma would compare
these terms without replacing the signed slope by its absolute value.

### 15.2 Analyze the signed work for \(r>2\)

For

\[
\varepsilon(t)=(1+t)^{-r},
\qquad r>2,
\]

the square-kick budget is finite. The remaining target is a block estimate for

\[
\sum_q\varepsilon(\log q)J_q\widehat B'(\log q^-)
-\int\varepsilon(t)e^{t/2}\widehat B'(t)\,dt.
\]

These are the discrete and continuous pieces of the same centered forcing.
An Abel-summation identity that keeps them together is more promising than
bounding them separately.

### 15.3 Use cell maxima, not dense grids

For \(\widehat B\), compute the exact Lambert-\(W\) point in each cell. For the
full \(B\), use a validated derivative-root solve. This makes any future
interval certification event-driven and avoids arbitrary sampling grids.

### 15.4 Certify only after an asymptotic lemma exists

The finite initial cone at \(q=2\) and any finite bridge can be certified by
directed interval arithmetic. Extending the current binary64 range alone does
not address infinity. The next major effort should therefore target the block
or kick inequality, not merely a larger sieve.

## 16. Relation to the literature

Masatoshi Suzuki's screw function provides a closely related RH-equivalent
positivity criterion. His function has the zero-side form

\[
\Psi_S(t)=\sum_\gamma\frac{1-e^{i\gamma t}}{\gamma^2},
\]

and Theorem 1.7 states that its pointwise non-negativity is equivalent to RH.
The present function

\[
\Phi(t)=\sum_\gamma
\frac{1-e^{i\gamma t}}{\gamma^2+1/4}
\]

uses a resolvent-shifted weight. Thus the one-sided positivity idea is part of
the known Weil/screw landscape; the specific prime-power cell ODE, event map,
and candidate cone above should be treated as a concrete reformulation and
research route, not as an established novelty claim.

Primary references:

- Masatoshi Suzuki, [Aspects of the screw function corresponding to the
  Riemann zeta function](https://arxiv.org/abs/2206.03682), especially
  Theorem 1.7 and the discussion of the Weil distribution.
- Masatoshi Suzuki, [The screw line of the Riemann zeta-function and its
  applications](https://arxiv.org/abs/2209.04658).
- Local derivation of \(B\):
  `working_notes/Rinf_Integral_Signed_Cancellation_Investigation_2026-08-04.md`.

## 17. Bottom line

The proposed branch was worth pursuing. It produced an exact prime-power
dynamical system and a substantially sharper proof target than a generic
request for tighter prime-counting envelopes.

The strongest current route is:

\[
\boxed{
\text{prove the prime-power kick inequality}
\quad
R(\log q)+\Phi'(\log q^-)
>\frac{\Lambda(q)}{\sqrt q}
\quad\text{for all }q=p^m.
}
\]

The full finite scan supports this invariant, and the open-cell part is already
exactly controlled. The missing universal control is confined to the discrete
prime-power kicks. That is real progress in localization, but it is still the
RH-strength step.
