from util import clamp_scalar, modinv

"""
A simple, from-scratch X25519 implementation in Python using only
built-in big integer support and the Montgomery ladder from RFC 7748.

WARNING: This code is not constant-time and is unsuitable for
production use. It is intended only as a pedagogical reference.
"""

# Curve25519 prime
P: int = 2**255 - 19

# A24 = (486662 - 2) / 4 for Curve25519
A24 = 121665


def decode_u(u_bytes: bytes) -> int:
    """
    Mask off the top bit to ensure compatibility with
    certain legacy encodings.
    """
    u = bytearray(u_bytes)
    u[31] &= 127
    return int.from_bytes(u, "little")


def encode_u(u_int: int) -> bytes:
    return u_int.to_bytes(32, "little")


def x25519_ladder(k_int: int, u_int: int) -> int:
    """
    Perform the Montgomery ladder (scalar multiplication) on Curve25519:
        X25519(k, u) = k * (u : 1) in the group law, returning the
        resulting u-coordinate as an integer.

    Follows the pseudo-code in RFC 7748, section 5.
    """
    # Montgomery ladder initialization
    x1 = u_int
    x2, z2 = 1, 0
    x3, z3 = u_int, 1
    swap = 0

    # Loop over bits of k from top (254) down to 0
    # (Note the top bit is always set after clamping, so we skip bit 255)

    # TODO: this works, but we can do better, e.g. look at Martin's tutorial

    for t in reversed(range(255)):
        k_t = (k_int >> t) & 1
        # swap ^= k_t
        # But more directly we can do a conditional swap:
        if k_t != swap:
            x2, x3 = x3, x2
            z2, z3 = z3, z2
            swap = k_t

        # The curve arithmetic
        #  -- all operations mod P
        A = (x2 + z2) % P
        B = (x2 - z2) % P
        AA = (A * A) % P
        BB = (B * B) % P
        E = (AA - BB) % P
        C = (x3 + z3) % P
        D = (x3 - z3) % P
        DA = (D * A) % P
        CB = (C * B) % P
        x3 = (DA + CB) % P
        x3 = (x3 * x3) % P
        z3 = (DA - CB) % P
        z3 = (z3 * z3) % P
        z3 = (z3 * x1) % P
        x2 = (AA * BB) % P
        z2 = (E * ((AA + (A24 * E) % P) % P)) % P

    # Last swap if needed
    if swap:
        x2, x3 = x3, x2
        z2, z3 = z3, z2

    # Return x2 * (z2^(p-2)) mod p  (the u-coordinate)
    # We use Fermat's little theorem for the inverse: z2^(p-2) mod p
    # inv_z2 = pow(z2, P - 2, P)
    inv_z2: int = modinv(z2, P)
    return (x2 * inv_z2) % P


def x25519(k_bytes: bytes, u_bytes: bytes) -> bytes:
    """
    Public, user-facing X25519 function:
      - Clamps scalar
      - Decodes u-coordinate
      - Montgomery ladder
      - Returns 32-byte little-endian result
    """
    k_int = clamp_scalar(bytearray(k_bytes))
    u_int = decode_u(u_bytes)
    x_result = x25519_ladder(k_int, u_int)
    return encode_u(x_result)


# --- TEST / DEMO CODE BELOW ---
if __name__ == "__main__":
    # Quick demonstration using the basepoint (9 + 31 zeros)

    # Example scalar from RFC 7748 test vectors
    k_hex = "a546e36bf0527c9d3b16154b82465edd" "62144c0ac1fc5a18506a2244ba449ac4"
    u_hex = "e6db6867583030db3594c1a424b15f7c" "726624ec26b3353b10a903a6d0ab1c4c"
    expect_hex = "c3da55379de9c6908e94ea4df28d084f" "32eccf03491c71f754b4075577a28552"

    k_bytes_ = bytes.fromhex(k_hex)
    u_bytes_ = bytes.fromhex(u_hex)
    result = x25519(k_bytes_, u_bytes_)

    # Demonstrate scalar multiplication with basepoint
    # (X25519 private key -> public key).
    basepoint = b"\x09" + b"\x00" * 31
    # Typical random 32-byte scalar:
    import os

    private_scalar = os.urandom(32)  # not cryptographically safe in example
    public_u = x25519(private_scalar, basepoint)

    public_u = x25519(private_scalar, basepoint)
