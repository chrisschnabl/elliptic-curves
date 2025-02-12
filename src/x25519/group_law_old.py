# Prime for Curve25519 and the curve parameter A.
from tonelli_shanks import tonelli
from util import clamp_scalar, encode_u_coordinate, modinv


# We represent the point at infinity as None.
Identity = None
CurvePoint = Point | None

# TODO: why is this different than the one in util.py?
def decode_public_key(pk: bytes) -> int:
    """
    Decodes a 32-byte X25519 public key into an integer u-coordinate.
    RFC 7748 requires that the most significant bit (bit 255) be cleared.
    """
    if len(pk) != 32:
        raise ValueError("Public key must be 32 bytes long")
    # Create a mutable copy and clear bit 255 of the last byte.
    pk_bytes = bytearray(pk)
    pk_bytes[31] &= 0x7F  # mask with 0x7F to clear the top bit
    # Convert the little-endian byte array to an integer.
    return int.from_bytes(pk_bytes, byteorder="little")

def recover_point(x: int) -> Point:
    """
    Given an x-coordinate, recover a point (x,y) on the curve defined by:
         y² = x³ + A*x² + x   (mod p)
    by computing a square root of the right-hand side.
    (If there is no square root, a ValueError is raised.)
    Note: In our variant the full point (x and y) is transmitted, so this function is not used.
    """
    rhs = (pow(x, 3, p) + A * pow(x, 2, p) + x) % p
    y = tonelli(rhs, p)
    if y is None:
        raise ValueError("No valid y for given x")
    # Choose the smaller square root.
    if y > p - y:
        y = p - y
    return (x, y)


def point_add(P: CurvePoint, Q: CurvePoint) -> CurvePoint:
    """
    Add two points P and Q (affine coordinates) on the Montgomery curve:
         y² = x³ + A*x² + x   (mod p)
    For distinct points:
         λ = (y₂ - y₁) / (x₂ - x₁) mod p,
         x₃ = λ² - A - x₁ - x₂ mod p,
         y₃ = λ*(x₁ - x₃) - y₁ mod p.
    If P == Q then doubling is performed.
    Returns the resulting point, or None if the result is the point at infinity.
    """
    if P is Identity:
        return Q
    if Q is Identity:
        return P

    x1, y1 = P
    x2, y2 = Q

    if x1 == x2:
        # If y₁ + y₂ ≡ 0 mod p, then Q is the inverse of P.
        if (y1 + y2) % p == 0:
            return Identity
        # Otherwise, P == Q and we perform doubling.
        return point_double(P)

    lam = ((y2 - y1) * modinv(x2 - x1, p)) % p
    x3 = (lam * lam - A - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)


def point_double(P: CurvePoint) -> CurvePoint:
    """
    Double the point P on the Montgomery curve.
    Uses the tangent line at P:
         λ = (3*x² + 2*A*x + 1) / (2*y) mod p,
         x₃ = λ² - A - 2*x mod p,
         y₃ = λ*(x - x₃) - y mod p.
    Returns the doubled point, or None if y == 0.
    """
    if P is Identity:
        return Identity
    x, y = P
    if y == 0:
        return Identity
    lam = ((3 * x * x + 2 * A * x + 1) * modinv(2 * y, p)) % p
    x3 = (lam * lam - A - 2 * x) % p
    y3 = (lam * (x - x3) - y) % p
    return (x3, y3)


def scalar_mult(k: int, P: CurvePoint) -> CurvePoint:
    """
    Compute the scalar multiplication k * P using a recursive double-and-add algorithm:
      - If k == 1, return P.
      - If k is even, return point_double(scalar_mult(k//2, P)).
      - If k is odd and k > 1, return point_add(point_double(scalar_mult((k-1)//2, P)), P).
    This method requires O(log k) doublings/additions.
    """
    if k == 1:
        return P
    if k % 2 == 0:
        return point_double(scalar_mult(k // 2, P))
    return point_add(point_double(scalar_mult((k - 1) // 2, P)), P)


def x25519(private_key_bytes: bytes, public_key_bytes: bytes) -> tuple[bytes, bytes]:
    """
    Compute the z25519 function using the group law on Curve25519.
    This function expects:
      - a 32-byte private key,
      - a 64-byte public key containing both the x-coordinate (u-coordinate) and the y-coordinate (each 32 bytes, little-endian).
    It returns a tuple (x_bytes, y_bytes) representing the resulting point.
    """
    scalar = clamp_scalar(bytearray(private_key_bytes))
    P = recover_point(decode_public_key(public_key_bytes))
    Q = scalar_mult(scalar, P)
    if Q is None:
        raise ValueError("Resulting point is the point at infinity")
    return encode_u_coordinate(Q[0])