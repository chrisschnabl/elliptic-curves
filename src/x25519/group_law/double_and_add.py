from util import clamp_scalar, modinv
from tonelli_shanks import tonelli

# Prime for Curve25519 and the curve parameter A.
p = 2**255 - 19
A = 486662
Point = tuple[int, int]

# Standard base point for Curve25519 (full coordinates).
BasePoint: Point = (
    9,
    #14781619447589544791020593568409986887264606134616475288964881837755586237401,
    43114425171068552920764898935933967039370386198203806730763910166200978582548,
)

# We represent the point at infinity as None.
MontgomeryIdentity = None
MontgomeryPoint = Point | None


def decode_public_key_old(pk: bytes) -> Point:
    """
    Decode a public key that includes both x and y coordinates.
    Expects a 64-byte input:
      - the first 32 bytes encode the x-coordinate (little-endian),
      - the next 32 bytes encode the y-coordinate (little-endian).
    Both coordinates are reduced modulo p.
    """
    #if len(pk) != 64:
    #    raise ValueError("Public key must be 64 bytes (32 for x, 32 for y)")
    x = int.from_bytes(pk[:32], "little") % p
    y = int.from_bytes(pk[32:], "little") % p
    return (x, y)

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

def decode_u(u_bytes: bytes) -> int:
    u = bytearray(u_bytes)
    u[31] &= 127
    return int.from_bytes(u, "little")

def encode_u_coordinate(x: int) -> bytes:
    """Encode an integer x as a 32-byte little-endian byte string."""
    return x.to_bytes(32, "little")


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


def point_add(P: MontgomeryPoint, Q: MontgomeryPoint) -> MontgomeryPoint:
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
    if P is MontgomeryIdentity:
        return Q
    if Q is MontgomeryIdentity:
        return P

    x1, y1 = P
    x2, y2 = Q

    if x1 == x2:
        # If y₁ + y₂ ≡ 0 mod p, then Q is the inverse of P.
        if (y1 + y2) % p == 0:
            return MontgomeryIdentity
        # Otherwise, P == Q and we perform doubling.
        return point_double(P)

    lam = ((y2 - y1) * modinv(x2 - x1, p)) % p
    x3 = (lam * lam - A - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)


def point_double(P: MontgomeryPoint) -> MontgomeryPoint:
    """
    Double the point P on the Montgomery curve.
    Uses the tangent line at P:
         λ = (3*x² + 2*A*x + 1) / (2*y) mod p,
         x₃ = λ² - A - 2*x mod p,
         y₃ = λ*(x - x₃) - y mod p.
    Returns the doubled point, or None if y == 0.
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

def main():
    # --- Test Vector 1 ---
    # The test vector provides a 32-byte scalar and a 32-byte u-coordinate.
    # Since in the group-law variant the full public key includes both x and y,
    # we recover y using recover_point.
    scalar_hex1    = "a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4"
    u_hex1         = "e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c"
    expected_u_hex1 = "c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552"

    scalar_bytes1 = bytes.fromhex(scalar_hex1)
    u_bytes1      = bytes.fromhex(u_hex1)
    # Recover the y-coordinate from the given x (u-coordinate).
    x1_int = int.from_bytes(u_bytes1, "little") % p
    # recover_point returns (x, y) with the smaller square root chosen.
    P1 = recover_point(x1_int)
    public_key_bytes1 = u_bytes1 + encode_u_coordinate(P1[1])

    result1 = x25519(scalar_bytes1, public_key_bytes1)

    print("Test Vector 1:")
    print("Input scalar (hex):         ", scalar_hex1)
    print("Input u-coordinate (hex):   ", u_hex1)
    print("Expected output u-coordinate (hex):", expected_u_hex1)
    print("Computed output u-coordinate (hex):", result1[0].hex())
    print()

    # --- Test Vector 2 ---
    scalar_hex2    = "4b66e9d4d1b4673c5ad22691957d6af5c11b6421e0ea01d42ca4169e7918ba0d"
    u_hex2         = "e5210f12786811d3f4b7959d0538ae2c31dbe7106fc03c3efc4cd549c715a493"
    expected_u_hex2 = "95cbde9476e8907d7aade45cb4b873f88b595a68799fa152e6f8f7647aac7957"

    scalar_bytes2 = bytes.fromhex(scalar_hex2)
    u_bytes2      = bytes.fromhex(u_hex2)
    x2_int = int.from_bytes(u_bytes2, "little") % p
    P2 = recover_point(x2_int)
    public_key_bytes2 = u_bytes2 + encode_u_coordinate(P2[1])

    result2 = z25519(scalar_bytes2, public_key_bytes2)

    print("Test Vector 2:")
    print("Input scalar (hex):         ", scalar_hex2)
    print("Input u-coordinate (hex):   ", u_hex2)
    print("Expected output u-coordinate (hex):", expected_u_hex2)
    print("Computed output u-coordinate (hex):", result2[0].hex())

if __name__ == "__main__":
    main()