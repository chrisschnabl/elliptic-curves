from util import clamp_scalar, modinv, tonelli

# Prime for Curve25519 and the curve parameter A.
p = 2**255 - 19
A = 486662
Point = tuple[int, int]

BasePoint: Point = (
    9,
    14781619447589544791020593568409986887264606134616475288964881837755586237401,
)


# TODO CS: here we can assume that the y is sent over the network
# TODO CS: so we don't need to recover it
# TODO the notes mention this but the slides ar unambigious


def recover_point(x: int) -> Point:
    """
    Given an x-coordinate (integer), recover an affine point (x,y) on the curve
    defined by y^2 = x^3 + A*x^2 + x (mod p).
    If no square root exists, raises ValueError.
    """
    rhs = (pow(x, 3, p) + A * pow(x, 2, p) + x) % p

    # if rhs % 8 == 5:
    #    y = modinv(rhs, p)
    # else:
    #    y = tonelli_shanks(rhs, p)
    # y = tonelli_shanks(rhs, p)
    y = tonelli(rhs, p)
    if y is None:
        raise ValueError("No valid y for given x")

    # (You may choose one of the two square roots; here we choose the smaller one)
    if y > p - y:  # TODO CS: why is this above needed
        y = p - y
    return (x, y)


# TODO CS: check this with the assignment again
# TODO CS: this uses projective coordinates
# TODO CS: use projective coordinate more often

MontgomeryIdentity = None
MontgomeryPoint = Point | MontgomeryIdentity


def point_add(P: MontgomeryPoint, Q: MontgomeryPoint) -> MontgomeryPoint:
    """
    Add two points P and Q (affine coordinates) on the Montgomery curve:
      y^2 = x^3 + A*x^2 + x  (mod p)
    For distinct P and Q:
       λ = (y2 - y1)/(x2 - x1) mod p,
       x3 = λ^2 - A - x1 - x2  mod p,
       y3 = λ*(x1 - x3) - y1   mod p. compared to slides this is simplified since we alreadyk now x3
       avoids recomputation, intuitively lambda * (x1 - x3) gives us the difference in y direction
    For doubling, see point_double.
    Returns the resulting point, or None if the result is the point at infinity.
    """
    if P is MontgomeryIdentity:
        return Q
    if Q is MontgomeryIdentity:
        return P

    x1, y1 = P
    x2, y2 = Q

    if x1 == x2:
        # If y1 == y2 then the result is the point at infinity.
        if (y1 + y2) % p == 0:
            return MontgomeryIdentity
        # Otherwise, treat it as doubling.
        return point_double(P)

    lam = ((y2 - y1) * modinv(x2 - x1, p)) % p
    x3 = (lam * lam - A - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)


def point_double(P: MontgomeryPoint) -> MontgomeryPoint:
    """
    Double the point P on the Montgomery curve.
    Uses the tangent line at P:
       λ = (3*x^2 + 2*A*x + 1)/(2*y) mod p,
       x3 = λ^2 - A - 2*x  mod p,
       y3 = λ*(x - x3) - y  mod p.
    Returns the doubled point, or None if y is 0.
    """
    if P is MontgomeryIdentity:
        return MontgomeryIdentity
    x, y = P
    if y == 0:
        return MontgomeryIdentity
    lam = ((3 * x * x + 2 * A * x + 1) * modinv(2 * y, p)) % p
    x3 = (lam * lam - A - 2 * x) % p
    y3 = (lam * (x - x3) - y) % p
    return (x3, y3)


def scalar_mult(k: int, P: MontgomeryPoint) -> MontgomeryPoint:
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
    # Just assume the key is 64 bytes and includes both x and y coordinates
    # P = (x_coord, int.from_bytes(public_key_bytes[32:], "little"))
    Q = scalar_mult(scalar, P)
    # Return the x-coordinate of the resulting point.
    if Q is None:
        raise ValueError("Resulting point is the point at infinity")
    return encode_u_coordinate(Q[0])


# Example usage
if __name__ == "__main__":
    """
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
    """
