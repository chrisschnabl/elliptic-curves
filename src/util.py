def legendre(a, p):
    return pow(a, (p - 1) // 2, p)


def legendre_symbol(a, p):
    ls = legendre(a, p)
    return -1 if ls == p - 1 else ls


def is_quadratic_residue(n, p):
    if n % p == 0 or p == 2:
        return True
    return legendre_symbol(n, p) == 1


def find_nonsquare(p):
    """
    Find and return a quadratic non-residue modulo p.
    (That is, find z such that legendre_symbol(z, p) == -1.)
    """
    for z in range(2, p):
        if legendre_symbol(z, p) == -1:
            return z
    raise ValueError(f"Could not find a quadratic non-residue modulo {p}")


def decompose(n):
    """
    Decompose n as q * 2^s with q odd.
    Returns a tuple (q, s).
    """
    s = 0
    while n % 2 == 0:
        s += 1
        n //= 2
    return n, s


def tonelli(n, p) -> None | int:
    """
    Solve for r in r^2 ≡ n (mod p).
    Returns a list [r, p-r] if a solution exists, or an empty list if no solution exists.
    """
    # Trivial cases:
    if p == 2:
        return n % p
    if n == 0:
        return 0
    if not is_quadratic_residue(n, p):
        return None

    # If p ≡ 3 (mod 4) then we can compute the solution directly.
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)

    # Factor p-1 as q * 2^s with q odd.
    q, s = decompose(p - 1)

    # Find a quadratic non-residue z modulo p.
    z = find_nonsquare(p)

    # Set up the initial values:
    c = pow(z, q, p)
    r = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s

    # Main loop: adjust r, t, c until t becomes 1.
    while t != 1:
        # Find the smallest i (0 < i < m) such that t^(2^i) ≡ 1 (mod p)
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


def tonelli_shanks(n: int, p: int) -> int | None:
    """
    Compute the square root of n modulo p (if it exists) using Tonelli-Shanks.
    Returns one square root, or None if no square root exists.
    """
    if n % p == 0:
        return 0
    # Check if n is a quadratic residue via Euler's criterion.
    if legendre(n, p) != 1:
        return None

    # Factor p - 1 as Q * 2^S with Q odd.
    S = 0
    Q = p - 1
    while Q % 2 == 0:
        Q //= 2
        S += 1

    if S == 1:
        return pow(n, (p + 1) // 4, p)

    # Find a quadratic non-residue z.
    z = 2
    while legendre(z, p) != p - 1:
        z += 1

    M = S
    c = pow(z, Q, p)
    t = pow(n, Q, p)
    R = pow(n, (Q + 1) // 2, p)

    while t != 1:
        # Find the smallest i (0 < i < M) such that t^(2^i) == 1 mod p.
        i = 0
        temp = t
        while temp != 1:
            temp = (temp * temp) % p
            i += 1
            if i == M:
                return None  # Should not happen
        b = pow(c, 2 ** (M - i - 1), p)
        M = i
        c = (b * b) % p
        t = (t * c) % p
        R = (R * b) % p
    return R


def clamp_scalar(k: bytearray) -> int:
    """
    Clamp the 32-byte scalar as per X25519 specification:
      - Clear the 3 least-significant bits.
      - Clear the most-significant bit.
      - Set the second-most-significant bit.

    Returns an integer (little-endian).
    """

    k[0] &= 248
    k[31] &= 127
    k[31] |= 64
    return int.from_bytes(k, "little")


def modinv(x: int, p: int) -> int:
    """Modular inverse modulo p (p is prime)."""
    return pow(x, p - 2, p)


def sqrt_mod(a: int, p: int) -> int:
    """
    Compute a square root of a modulo p using the following method:
      Let r = a^((p+3)//8) mod p.
      If r^2 ≡ a mod p then return r.
      Else if r^2 ≡ -a mod p then return (r * sqrt(-1)) mod p,
         where sqrt(-1) = 2^((p-1)//4) mod p.
      Otherwise, raise an error.

    This method works for p = 2^255 - 19 (I guess)
    """

    exp = (p + 3) // 8
    r = pow(a, exp, p)
    if (r * r) % p == a % p:
        return r
    if (r * r) % p == (-a) % p:
        sqrt_m1 = pow(2, (p - 1) // 4, p)
        return (r * sqrt_m1) % p
    raise ValueError("No square root exists for the given input.")


def decode_public_key(u_bytes: bytes, p: int) -> int:
    """
    Decode a 32-byte public key: interpret as little-endian integer,
    clear the most significant bit, then reduce modulo p.
    """
    u_int = int.from_bytes(u_bytes, "little")
    u_int &= (1 << 255) - 1
    return u_int % p


def encode_u_coordinate(x: int) -> bytes:
    """Encode an integer x as a 32-byte little-endian byte string."""
    return x.to_bytes(32, "little")
