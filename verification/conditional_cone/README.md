# Finite verification for the strict prime-power cone

The accompanying paper uses a directed interval calculation on $2\le x\le60\,000$ and an RH-conditional analytic argument on $x\ge e^{11}$. This script reproduces the finite calculation without any assumption about the zeros of the Riemann zeta function.

Use Python 3.10 or later. From the root of the research archive, run

    python -B verification/conditional_cone/verify_cone.py --check

The command recomputes the intervals, compares the output with `report.json`, and exits with a nonzero status on a mismatch or a failed positivity check. Running it without options prints the recomputed report. The only project dependency is [the August Decimal interval backend](../finite_robin_prime_power/test_ca_residual_finite_certificate.py). Keep that file in its original sibling directory when copying the computation elsewhere. Both files use only the Python standard library.

The sieve marks the base $p$ of every integer prime power $q=p^k\le60\,000$. At each such $q$, in increasing order, the script adds $\log p$ to $\psi(q)$ and $(\log p)/q$ to $S(q)$. Both sums include the current event, fixing the right-derivative convention. With $C=2+\gamma-\log(4\pi)$ it evaluates

$$
\begin{aligned}
u &= S(q)-\log q+\gamma+q^{-1}-\operatorname{atanh}(q^{-1}),\\
v &= \psi(q)-q+\log(2\pi)+\tfrac12\log(1-q^{-2}),\\
\Phi &= C-\sqrt q\,u+v/\sqrt q,\qquad
W=\tfrac12(\sqrt q\,u+v/\sqrt q),\\
k &= \sqrt q-\frac1{\sqrt q(q^2-1)}-\frac C4,\qquad
D=\Phi(2k+\Phi/4)-W^2.
\end{aligned}
$$

The inequalities checked are $\inf\Phi>0$ and $\inf D>0$ immediately after every event. The smallest lower endpoints in the archived run exceed $0.0277660630$ for $\Phi$, at $q=3109$, and $0.0789828268$ for $D$, at $q=2$. All 6,148 events are covered. The paper extends the endpoint enclosures across the open cells, including the final segment up to 60,000, using $D'=2k'\Phi$. The script also checks that the upper endpoint of the interval for $e^{11}$ lies below 60,000.

Basic operations use floor and ceiling Decimal contexts with 50 significant digits. The backend encloses logarithms and square roots by padding correctly rounded results outwards, $\pi$ by Machin's formula and $\gamma$ by harmonic-sum inequalities with 100,000 terms. The certificate relies on these enclosures, the sieve and the stated equations of the paper.
