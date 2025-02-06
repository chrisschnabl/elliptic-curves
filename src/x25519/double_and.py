import os

from x25519.util import clamp_scalar, modinv

# Prime for Curve25519 and the curve parameter A.
p = 2**255 - 19
A = 486662
Point = tuple[int, int]


def tonelli_shanks(n: int, p: int) -> int | None:
    """
    Compute the square root of n modulo p (if it exists) using Tonelli-Shanks.
    Returns one square root, or None if no square root exists.
    """
    # Check if n is a quadratic residue via Euler's criterion.
    if pow(n, (p - 1) // 2, p) != 1:
        return None

    # Factor p - 1 as Q * 2^S with Q odd.
    S = 0
    Q = p - 1
    while Q % 2 == 0:
        Q //= 2
        S += 1

    # Find a quadratic non-residue z.
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
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


def recover_point(x: int) -> Point:
    """
    Given an x-coordinate (integer), recover an affine point (x,y) on the curve
    defined by y^2 = x^3 + A*x^2 + x (mod p).
    If no square root exists, raises ValueError.
    """
    rhs = (pow(x, 3, p) + A * pow(x, 2, p) + x) % p
    y = tonelli_shanks(rhs, p)
    if y is None:
        raise ValueError("No valid y for given x")
    # (You may choose one of the two square roots; here we choose the smaller one)
    if y > p - y:  # TODO CS: why is this above needed
        y = p - y
    return (x, y)


#  TODO CS: check this with the assignment again
# TODO CS: this uses projective coordinates
# TODO CS: use projective coordinate more often
def point_add(P: Point | None, Q: Point | None) -> Point | None:
    """
    Add two points P and Q (affine coordinates) on the Montgomery curve:
      y^2 = x^3 + A*x^2 + x  (mod p)
    For distinct P and Q:
       λ = (y2 - y1)/(x2 - x1) mod p,
       x3 = λ^2 - A - x1 - x2  mod p,
       y3 = λ*(x1 - x3) - y1   mod p.
    For doubling, see point_double.
    Returns the resulting point, or None if the result is the point at infinity.
    """
    if P is None:
        return Q
    if Q is None:
        return P

    x1, y1 = P
    x2, y2 = Q

    # Handle the case where the x-coordinates are equal.
    if x1 == x2:
        # If y1 is the negative of y2 then the result is the point at infinity.
        if (y1 + y2) % p == 0:
            return None
        # Otherwise, treat it as doubling.
        return point_double(P)

    lam = ((y2 - y1) * modinv(x2 - x1, p)) % p
    x3 = (lam * lam - A - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)


def point_double(P: Point | None) -> Point | None:
    """
    Double the point P on the Montgomery curve.
    Uses the tangent line at P:
       λ = (3*x^2 + 2*A*x + 1)/(2*y) mod p,
       x3 = λ^2 - A - 2*x  mod p,
       y3 = λ*(x - x3) - y  mod p.
    Returns the doubled point, or None if y is 0.
    """
    if P is None:
        return None
    x, y = P
    if y == 0:
        return None
    lam = ((3 * x * x + 2 * A * x + 1) * modinv(2 * y, p)) % p
    x3 = (lam * lam - A - 2 * x) % p
    y3 = (lam * (x - x3) - y) % p
    return (x3, y3)


def scalar_mult(k: int, P: Point | None) -> Point | None:
    """
    Compute the scalar multiplication k * P using the recursive double-and-add algorithm:
      - If k == 1, return P.
      - If k is even, return point_double(scalar_mult(k//2, P)).
      - If k is odd and k > 1, ret point_add(point_double(scalar_mult((k-1)//2, P)), P).
    This computes k*P in O(log k) additions/doublings.
    """
    if k == 1:
        return P
    if k % 2 == 0:
        return point_double(scalar_mult(k // 2, P))
    return point_add(point_double(scalar_mult((k - 1) // 2, P)), P)


def decode_public_key(u_bytes: bytes) -> int:
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


def x25519(private_key_bytes: bytes, public_key_bytes: bytes) -> bytes:
    """
    Compute the X25519 function:
      - Clamp the 32-byte private key.
      - Decode the public key (interpreted as the x-coordinate of a point).
      - Recover the corresponding affine point P on Curve25519.
      - Compute Q = (clamped scalar) * P using double-and-add scalar multiplication.
      - Return the x-coordinate of Q as 32 bytes (little-endian).
    """
    scalar = clamp_scalar(bytearray(private_key_bytes))
    x_coord = decode_public_key(public_key_bytes)
    P = recover_point(x_coord)
    # TODO CS: why do we need this, and can we make this easier
    # TODO CS: why do we have to work in affine points
    Q = scalar_mult(scalar, P)
    # Return the x-coordinate of the resulting point.
    if Q is None:
        raise ValueError("Resulting point is the point at infinity")
    return encode_u_coordinate(Q[0])


# Example usage
if __name__ == "__main__":
    # Generate a random 32-byte private key.
    private_key = os.urandom(32)
    # For demonstration, use the standard base point whose x-coordinate is 9.
    base_x = 9
    public_key = base_x.to_bytes(32, "little")  # This is our input “public key”

    # Compute "public key" = private_key * base_point, using our function.
    result = x25519(private_key, public_key)

    # Quick demonstration using the basepoint (9 + 31 zeros)

    # Example scalar from RFC 7748 test vectors
    k_hex = "a546e36bf0527c9d3b16154b82465edd" "62144c0ac1fc5a18506a2244ba449ac4"
    u_hex = "e6db6867583030db3594c1a424b15f7c" "726624ec26b3353b10a903a6d0ab1c4c"
    expect_hex = "c3da55379de9c6908e94ea4df28d084f" "32eccf03491c71f754b4075577a28552"

    k_bytes_ = bytes.fromhex(k_hex)
    u_bytes_ = bytes.fromhex(u_hex)
    result = x25519(k_bytes_, u_bytes_)
