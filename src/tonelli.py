from typing import Optional, Tuple

# Translated into Python from postscript code found here:
# https://github.com/mvaneerde/blog/blob/develop/tonelli-shanks/tonelli-shanks.ps1

def decompose(n: int) -> Tuple[int, int]:
    """
    Decompose an integer n as q * 2^s with q odd.

    Args:
        n: A positive integer.

    Returns:
        A tuple (q, s) such that n = q * 2^s and q is odd.
    """
    s = 0
    while n % 2 == 0:
        s += 1
        n //= 2
    return n, s  # Here, n is the odd factor q.


def legendre(a: int, p: int) -> int:
    """
    Compute the raw Legendre symbol value (as a power modulo p).

    Args:
        a: The integer whose Legendre symbol is to be computed.
        p: An odd prime.

    Returns:
        pow(a, (p - 1) // 2, p) which is either 1, 0, or p-1.
    """
    return pow(a, (p - 1) // 2, p)


def legendre_symbol(a: int, p: int) -> int:
    """
    Compute the Legendre symbol (a/p).

    Args:
        a: The integer to test.
        p: An odd prime.

    Returns:
        1 if a is a quadratic residue modulo p,
       -1 if a is a non-residue,
        0 if a is divisible by p.
    """
    ls = legendre(a, p)
    return -1 if ls == p - 1 else ls


def is_quadratic_residue(n: int, p: int) -> bool:
    """
    Check if n is a quadratic residue modulo p.

    Args:
        n: The integer to test.
        p: A prime number.

    Returns:
        True if n is a quadratic residue modulo p, or if n is divisible by p
        or if p == 2; otherwise, False.
    """
    if n % p == 0 or p == 2:
        return True
    return legendre_symbol(n, p) == 1


def find_nonsquare(p: int) -> int:
    """
    Find a quadratic non-residue modulo an odd prime p.

    That is, find z in the range 2 <= z < p such that legendre_symbol(z, p) == -1.

    Args:
        p: An odd prime.

    Returns:
        An integer z that is a quadratic non-residue modulo p.

    Raises:
        ValueError: If no quadratic non-residue is found (should not occur for prime p).
    """
    for z in range(2, p):
        if legendre_symbol(z, p) == -1:
            return z
    raise ValueError(f"Could not find a quadratic non-residue modulo {p}")


def tonelli(n: int, p: int) -> Optional[int]:
    """
    Solve for a square root r of n modulo p, i.e. find r such that r^2 ≡ n (mod p).

    For an odd prime p and a quadratic residue n modulo p, there are two solutions: r and p - r.
    This function returns one of the solutions. (The other solution can be obtained as p - r.)
    For p == 2, the solution is trivially n mod 2.

    Args:
        n: The integer whose square root modulo p is to be computed.
        p: An odd prime (or p == 2).

    Returns:
        A square root r (0 <= r < p) if one exists, or None if n is not a quadratic residue modulo p.
    """
    # Special cases
    if p == 2:
        return n % 2
    if n == 0:
        return 0
    if not is_quadratic_residue(n, p):
        return None

    # For primes where p % 4 == 3, a direct solution exists.
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)

    # Factor p-1 as q * 2^s with q odd.
    q, s = decompose(p - 1)

    # Find a quadratic non-residue z modulo p.
    z = find_nonsquare(p)

    # Initialize variables for the Tonelli-Shanks algorithm.
    c = pow(z, q, p)
    r = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s

    # Main loop: update r, t, c until t becomes 1.
    while t != 1:
        # Find the smallest integer i (0 < i < m) such that t^(2^i) ≡ 1 (mod p).
        i = 0
        temp = t
        while temp != 1:
            temp = pow(temp, 2, p)
            i += 1
            if i == m:
                raise ValueError("Algorithm error: t^(2^i) never reached 1")

        # Compute b = c^(2^(m-i-1)) mod p.
        b = pow(c, 1 << (m - i - 1), p)

        # Update r, t, c, and m.
        r = (r * b) % p
        t = (t * b * b) % p
        c = (b * b) % p
        m = i

    return r