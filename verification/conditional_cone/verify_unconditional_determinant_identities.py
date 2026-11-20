"""Exact rational checks for the determinant and integer-factorial identities."""

from collections import Counter
from fractions import Fraction
from itertools import product
from math import factorial, isqrt


def check_polynomials():
    # After clearing X^2=q, each variable has degree at most two.
    # Three distinct values per variable suffice for polynomial equality.
    grid = tuple(map(Fraction, (-1, 0, 1)))
    for xroot, k, origin, z, v in product(grid, repeat=5):
        # z = F - log(2*pi) - q*T(1/q)/2 = X*B; X need not be divided out.
        xp = origin * xroot - z
        xw = z/2 + v
        scaled_d = xp * (2*k*xroot + xp/4) - xw*xw
        coefficient = (2*k + origin/2)*xroot + v
        if scaled_d != coefficient*xp - (v + origin*xroot/2)**2:
            raise ArithmeticError("The completed square has the wrong coefficients")

    def determinant(phi, k, slope):
        return phi * (2*k + phi/4) - slope**2

    for phi, k, slope, shift in product(grid, repeat=4):
        delta = determinant(phi+shift, k-shift/4, slope) - determinant(phi, k, slope)
        if delta != 2*shift*k - shift**2/4:
            raise ArithmeticError("Incorrect shifted determinant")


def factor(n):
    """Factor the next integer independently of the prime-power summation."""
    result = Counter()
    for p in range(2, isqrt(n)+1):
        while n % p == 0:
            result[p] += 1
            n //= p
    if n > 1:
        result[n] += 1
    return result


def check_factorial(limit):
    # Keep rational coefficients of each log(p), without evaluating the logs.
    primes = [p for p in range(2, limit+1) if all(p % d for d in range(2, isqrt(p)+1))]
    factorial_valuations = Counter()
    previous_a, previous_s = Counter(), Counter()
    for n in range(2, limit+1):
        increments = factor(n)
        factorial_valuations.update(increments)
        a, reciprocal, ramp = Counter(), Counter(), Counter()
        for p in primes:
            if p > n:
                break
            power = p
            while power <= n:
                a[p] += 1 - Fraction(n % power, power)
                reciprocal[p] += Fraction(1, power)
                ramp[p] += Fraction(n, power) - 1
                power *= p
            if ramp[p] + a[p] != factorial_valuations[p]:
                raise ArithmeticError(f"Factorial decomposition failed at N={n}, p={p}")
            if a[p]-previous_a[p] != increments[p]-previous_s[p]:
                raise ArithmeticError(f"First difference failed at N={n}, p={p}")
        previous_a, previous_s = a, reciprocal


def check_enlarged_factorial(q_limit, m_limit):
    """Compare coefficients of M*A(q) in the four-term enlarged-factorial identity."""
    maximum = q_limit * m_limit
    primes = [p for p in range(2, maximum+1) if all(p % d for d in range(2, isqrt(p)+1))]
    powers = {}
    for p in primes:
        powers[p] = []
        n = p
        while n <= maximum:
            powers[p].append(n)
            n *= p

    def factorial_exponents(n):
        # Factor the actual integer n! without using Legendre's floor-sum formula.
        value = factorial(n)
        exponents = Counter()
        for p in primes:
            if p > n:
                break
            while value % p == 0:
                exponents[p] += 1
                value //= p
        if value != 1:
            raise ArithmeticError(f"Incomplete factorization of {n}!")
        return exponents

    cases = list(product(range(2, q_limit+1), range(2, m_limit+1)))
    needed = set(range(2, q_limit+1)) | {q*m for q, m in cases}
    valuations = {n: factorial_exponents(n) for n in needed}
    comparisons = 0
    for q, m in cases:
        mq = m*q
        for p in primes:
            if p > mq:
                break
            small = [n for n in powers[p] if n <= q]
            left = m*sum(1 - Fraction(q % n, n) for n in small)
            # For each psi(Mq/j), include all eligible powers, including the endpoint.
            psi_count = sum(j*n <= mq for j in range(1, m+1) for n in powers[p])
            fractional_tail = sum(Fraction(mq % n, n) for n in small)
            right = m*valuations[q][p] - valuations[mq][p] + psi_count - fractional_tail
            if left != right:
                raise ArithmeticError(f"Enlarged factorial failed at q={q}, M={m}, p={p}")
            comparisons += 1
    return len(cases), comparisons


if __name__ == "__main__":
    check_polynomials()
    print("PASS: cleared determinant identity (243 grid points); shift identity (81 grid points)")
    check_factorial(300)
    print("PASS: exact coefficients of log(p) in the factorial identity and recurrence, N=2..300")
    cases, comparisons = check_enlarged_factorial(64, 10)
    print(f"PASS: enlarged factorial identity, q=2..64, M=2..10: "
          f"{cases} pairs and {comparisons} exact coefficient comparisons")
