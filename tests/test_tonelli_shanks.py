import unittest

from util import tonelli
from x25519.double_and import tonelli_shanks


def primes_method5(n):
    out = []
    sieve = [True] * (n + 1)
    for p in range(2, n + 1):
        if sieve[p]:
            out.append(p)
            for i in range(p, n + 1, p):
                sieve[i] = False
    return out


def legendre(a, p):
    return pow(a, (p - 1) // 2, p)


def tonelli2(n, p):
    """
    Solve for r in Z/pZ with r^2 ≡ n (mod p),
    assuming that n is a quadratic residue modulo the odd prime p.

    Returns:
        A square root r with 0 ≤ r < p, or
        None if no square root exists.
    """

    # Special case: if n is 0, then the square root is 0.
    n %= p
    if n == 0:
        return 0

    # Check that n is a quadratic residue via the Legendre symbol.
    # (Assumes legendre(n, p) returns 1 if n is a residue and -1 if not.
    #  If your legendre returns p-1 in place of -1, change the test accordingly.)
    if legendre(n, p) != 1:
        return None  # n is not a quadratic residue modulo p.

    # Write p-1 as Q * 2^S with Q odd.
    Q = p - 1
    S = 0
    while Q % 2 == 0:
        Q //= 2
        S += 1

    # If S is 1, we can solve the congruence directly.
    if S == 1:
        return pow(n, (p + 1) // 4, p)

    # Find a quadratic non-residue z (i.e. legendre(z, p) == -1).
    z = 2
    while legendre(z, p) != -1:
        z += 1

    # Set initial values:
    # c = z^Q mod p, R = n^{(Q+1)/2} mod p, t = n^Q mod p, and M = S.
    c = pow(z, Q, p)
    R = pow(n, (Q + 1) // 2, p)
    t = pow(n, Q, p)
    M = S

    # Loop until t becomes 1.
    while t != 1:
        # Find the smallest integer i, 0 < i < M, such that t^(2^i) ≡ 1 (mod p).
        i = 0
        temp = t
        while temp != 1:
            temp = (temp * temp) % p
            i += 1
            if i == M:  # (should not happen if n is a residue)
                break

        # Compute b = c^(2^(M-i-1)) mod p.
        b = pow(c, 1 << (M - i - 1), p)
        # Update R, t, c and M.
        R = (R * b) % p
        t = (t * b * b) % p
        c = (b * b) % p
        M = i

    return R


class TestTonelliShanks(unittest.TestCase):
    def test_tonelli_shanks_against_naive(self):
        for n in range(2, 1000):
            p = 2**255 - 19
            ts = tonelli(n, p)
            t = ts**2 % p if ts is not None else None
            t2 = (
                tonelli_shanks(n, p) ** 2 % p
                if tonelli_shanks(n, p) is not None
                else None
            )
            self.assertEqual(t, t2, f"n: {n}, p: {p}")
