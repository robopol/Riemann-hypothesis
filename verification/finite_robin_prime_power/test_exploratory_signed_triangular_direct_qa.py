#!/usr/bin/env python3
"""Independent high-precision QA for the signed triangular identity.

This test deliberately does not import the exploratory scanner.  It rebuilds
the small prime-power support, evaluates Phi from finite prime-side sums, and
checks T=S-I against 2*h*(p_bar-D_h) using event-split Gauss-Legendre
quadrature.
"""

from __future__ import annotations

import bisect
import math
import unittest

import numpy as np

try:
    import mpmath as mp
except ImportError:  # pragma: no cover - optional high-precision dependency
    mp = None


@unittest.skipUnless(mp is not None, "mpmath is required for this QA test")
class SignedTriangularDirectQATests(unittest.TestCase):
    """Verify the signs and midpoint convention without scanner identities."""

    @staticmethod
    def primes_upto(limit: int) -> list[int]:
        """Return primes through limit using a standalone dense sieve."""

        flags = bytearray(b"\x01") * (limit + 1)
        flags[:2] = b"\x00\x00"
        for prime in range(2, math.isqrt(limit) + 1):
            if flags[prime]:
                start = prime * prime
                flags[start : limit + 1 : prime] = b"\x00" * (
                    ((limit - start) // prime) + 1
                )
        return [value for value, is_prime in enumerate(flags) if is_prime]

    @classmethod
    def event_support(cls, limit: int):
        """Build all q=p^a events independently of the production helper."""

        rows = []
        for prime in cls.primes_upto(limit):
            value = prime
            exponent = 1
            log_prime = mp.log(prime)
            while value <= limit:
                rows.append((value, prime, exponent, log_prime))
                value *= prime
                exponent += 1
        rows.sort(key=lambda row: row[0])
        return rows

    def test_direct_signed_triangular_identity_at_q_997(self) -> None:
        """Check T=S-I=2*h*E at q=997 and lambda=1/2."""

        mp.mp.dps = 80
        q = 997
        central_prime = 997
        lambda_value = mp.mpf("0.5")
        t = mp.log(q)
        h = lambda_value * t / mp.sqrt(q)
        support_limit = int(mp.ceil(q * mp.exp(h))) + 5
        events = self.event_support(support_limit)
        event_values = [row[0] for row in events]

        cumulative_psi = []
        cumulative_s_lambda = []
        psi = mp.mpf(0)
        s_lambda = mp.mpf(0)
        for value, _prime, _exponent, log_prime in events:
            psi += log_prime
            s_lambda += log_prime / value
            cumulative_psi.append(psi)
            cumulative_s_lambda.append(s_lambda)

        weil_c = 2 + mp.euler - mp.log(4 * mp.pi)
        log_two_pi = mp.log(2 * mp.pi)

        def state_at_x(x):
            position = bisect.bisect_right(event_values, x) - 1
            psi_x = cumulative_psi[position] if position >= 0 else mp.mpf(0)
            s_lambda_x = (
                cumulative_s_lambda[position] if position >= 0 else mp.mpf(0)
            )
            log_x = mp.log(x)
            root_x = mp.sqrt(x)
            y = 1 / x
            trivial_f = (
                (1 + y) * mp.log1p(y)
                + (y - 1) * mp.log1p(-y)
                - 2 * y
            )
            trivial_f_prime = mp.log1p(-(y * y))
            u = (
                s_lambda_x
                - log_x
                + mp.euler
                + y * trivial_f_prime / 2
                - trivial_f / 2
            )
            v = psi_x - x + log_two_pi + trivial_f_prime / 2
            phi = weil_c - (root_x * u - v / root_x)
            post_slope = -(root_x * u + v / root_x) / 2
            return phi, post_slope

        phi, post_slope = state_at_x(mp.mpf(q))
        jump = mp.log(central_prime) / mp.sqrt(q)
        midpoint_slope = post_slope + jump / 2
        phi_minus, _ = state_at_x(mp.exp(t - h))
        phi_plus, _ = state_at_x(mp.exp(t + h))
        central_secant = (phi_plus - phi_minus) / (2 * h)
        discrepancy = midpoint_slope - central_secant

        signed_atomic_sum = mp.mpf(0)
        breakpoints = [mp.mpf(0), h]
        for value, _prime, _exponent, log_prime in events:
            offset = mp.log(value) - t
            event_jump = log_prime / mp.sqrt(value)
            if 0 < offset < h:
                signed_atomic_sum += (h - offset) * event_jump
                breakpoints.append(offset)
            elif -h < offset < 0:
                signed_atomic_sum -= (h + offset) * event_jump
                breakpoints.append(-offset)

        def forcing(time_value):
            x = mp.exp(time_value)
            phi_value, _ = state_at_x(x)
            q_factor = 2 * mp.sqrt(x) * (1 - 1 / (x * (x * x - 1)))
            return phi_value / 4 + q_factor / 2 - weil_c / 4

        def regular_integral(order: int):
            nodes, weights = mp.gauss_quadrature(order, "legendre")
            total = mp.mpf(0)
            ordered = sorted(set(breakpoints))
            for left, right in zip(ordered, ordered[1:]):
                midpoint = (left + right) / 2
                half_width = (right - left) / 2
                for node, weight in zip(nodes, weights):
                    s = midpoint + half_width * node
                    integrand = (h - s) * (
                        forcing(t + s) - forcing(t - s)
                    )
                    total += half_width * weight * integrand
            return total

        integral_16 = regular_integral(16)
        integral_32 = regular_integral(32)
        triangular_term = signed_atomic_sum - integral_32
        identity_target = 2 * h * discrepancy

        self.assertLess(abs(integral_16 - integral_32), mp.mpf("1e-60"))
        self.assertLess(abs(triangular_term - identity_target), mp.mpf("1e-60"))
        self.assertAlmostEqual(
            float(signed_atomic_sum), 0.019867885748400842, places=16
        )
        self.assertAlmostEqual(float(integral_32), 0.006886950883731535, places=16)
        self.assertAlmostEqual(
            float(triangular_term), 0.012980934864669307, places=16
        )


class SignedTriangularDirectBinary64QATests(unittest.TestCase):
    """Run the independent event-split check without optional packages."""

    @staticmethod
    def primes_upto(limit: int) -> list[int]:
        """Return primes through limit using a standalone dense sieve."""

        flags = bytearray(b"\x01") * (limit + 1)
        flags[:2] = b"\x00\x00"
        for prime in range(2, math.isqrt(limit) + 1):
            if flags[prime]:
                start = prime * prime
                flags[start : limit + 1 : prime] = b"\x00" * (
                    ((limit - start) // prime) + 1
                )
        return [value for value, is_prime in enumerate(flags) if is_prime]

    def test_direct_identity_binary64_fallback(self) -> None:
        """Recompute T=S-I at q=997 with split Gauss-Legendre quadrature."""

        q = 997
        lambda_value = 0.5
        t = math.log(q)
        h = lambda_value * t / math.sqrt(q)
        support_limit = math.ceil(q * math.exp(h)) + 5
        events: list[tuple[int, int, float]] = []
        for prime in self.primes_upto(support_limit):
            value = prime
            exponent = 1
            log_prime = math.log(prime)
            while value <= support_limit:
                events.append((value, exponent, log_prime))
                value *= prime
                exponent += 1
        events.sort(key=lambda row: row[0])
        event_values = [row[0] for row in events]

        cumulative_psi: list[float] = []
        cumulative_s_lambda: list[float] = []
        psi = 0.0
        s_lambda = 0.0
        for value, _exponent, log_prime in events:
            psi += log_prime
            s_lambda += log_prime / value
            cumulative_psi.append(psi)
            cumulative_s_lambda.append(s_lambda)

        euler_gamma = 0.577215664901532860606512090082402431
        log_two_pi = math.log(2.0 * math.pi)
        weil_c = 2.0 + euler_gamma - math.log(4.0 * math.pi)

        def state_at_x(x: float) -> tuple[float, float]:
            position = bisect.bisect_right(event_values, x) - 1
            psi_x = cumulative_psi[position] if position >= 0 else 0.0
            s_lambda_x = (
                cumulative_s_lambda[position] if position >= 0 else 0.0
            )
            log_x = math.log(x)
            root_x = math.sqrt(x)
            y = 1.0 / x
            trivial_f = (
                (1.0 + y) * math.log1p(y)
                + (y - 1.0) * math.log1p(-y)
                - 2.0 * y
            )
            trivial_f_prime = math.log1p(-(y * y))
            u = (
                s_lambda_x
                - log_x
                + euler_gamma
                + y * trivial_f_prime / 2.0
                - trivial_f / 2.0
            )
            v = psi_x - x + log_two_pi + trivial_f_prime / 2.0
            phi = weil_c - (root_x * u - v / root_x)
            post_slope = -(root_x * u + v / root_x) / 2.0
            return phi, post_slope

        _phi, post_slope = state_at_x(float(q))
        jump = math.log(q) / math.sqrt(q)
        midpoint_slope = post_slope + jump / 2.0
        phi_minus, _ = state_at_x(math.exp(t - h))
        phi_plus, _ = state_at_x(math.exp(t + h))
        central_secant = (phi_plus - phi_minus) / (2.0 * h)
        discrepancy = midpoint_slope - central_secant

        signed_atomic_sum = 0.0
        breakpoints = [0.0, h]
        for value, _exponent, log_prime in events:
            offset = math.log(value) - t
            event_jump = log_prime / math.sqrt(value)
            if 0.0 < offset < h:
                signed_atomic_sum += (h - offset) * event_jump
                breakpoints.append(offset)
            elif -h < offset < 0.0:
                signed_atomic_sum -= (h + offset) * event_jump
                breakpoints.append(-offset)

        def forcing(time_value: float) -> float:
            x = math.exp(time_value)
            phi_value, _ = state_at_x(x)
            q_factor = 2.0 * math.sqrt(x) * (
                1.0 - 1.0 / (x * (x * x - 1.0))
            )
            return phi_value / 4.0 + q_factor / 2.0 - weil_c / 4.0

        nodes, weights = np.polynomial.legendre.leggauss(32)
        regular_integral = 0.0
        ordered = sorted(set(breakpoints))
        for left, right in zip(ordered, ordered[1:]):
            midpoint = (left + right) / 2.0
            half_width = (right - left) / 2.0
            for node, weight in zip(nodes, weights):
                s = midpoint + half_width * float(node)
                regular_integral += (
                    half_width
                    * float(weight)
                    * (h - s)
                    * (forcing(t + s) - forcing(t - s))
                )

        triangular_term = signed_atomic_sum - regular_integral
        identity_target = 2.0 * h * discrepancy
        self.assertLess(abs(triangular_term - identity_target), 5e-13)
        self.assertAlmostEqual(signed_atomic_sum, 0.019867885748400842, places=14)
        self.assertAlmostEqual(regular_integral, 0.006886950883731535, places=14)
        self.assertAlmostEqual(triangular_term, 0.012980934864669307, places=14)


if __name__ == "__main__":
    unittest.main()
