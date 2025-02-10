from util import clamp_scalar, modinv  # Assumed to be available

# --- Global constants for Curve25519 ---
P: int = 2**255 - 19
A24: int = 121665  # (486662 - 2) / 4

# --- Conversion Helpers ---

def affine_to_projective(x: int) -> (int, int):
    """
    Convert an affine coordinate to a projective coordinate.
    Given x (affine), return (X:Z) with Z = 1.
    """
    return (x, 1)

def projective_to_affine(X: int, Z: int) -> int:
    """
    Convert a projective coordinate (X:Z) to an affine coordinate x = X/Z mod P.
    """
    inv_Z = modinv(Z, P)
    return (X * inv_Z) % P

# --- Core Arithmetic Helpers ---

def cswap(swap: int, x: int, y: int) -> (int, int):
    """
    Branchless conditional swap.

    Given a flag (swap = 0 or 1) and two numbers x and y,
    return (x, y) unchanged if swap == 0 or swapped if swap == 1.
    (Note: in Python this is not constant‑time, but it shows the intended arithmetic.)
    """
    dummy = swap * (x - y)
    x_new = x - dummy
    y_new = y + dummy
    return x_new, y_new

def ladder_step(X0: int, Z0: int, X1: int, Z1: int, xP: int) -> (int, int, int, int):
    """
    Perform one simultaneous step of the Montgomery ladder:
    
      - R0 = (X0:Z0) is doubled using projective doubling.
      - R1 = (X1:Z1) is updated via differential addition with R0.
    
    The formulas are:
      * Doubling:
          A  = X0 + Z0
          B  = X0 - Z0
          AA = A^2,   BB = B^2
          C  = AA - BB
          new_X0 = AA * BB
          new_Z0 = C * (AA + A24 * C)
      
      * Differential addition (using the fact that R1 - R0 = (xP:1)):
          D = X1 + Z1
          E = X1 - Z1
          DA = E * A
          CB = D * B
          new_X1 = (DA + CB)^2
          new_Z1 = xP * (DA - CB)^2
          
    All operations are performed modulo P.
    """
    # --- Doubling of R0 ---
    A = (X0 + Z0) % P
    B = (X0 - Z0) % P
    AA = (A * A) % P
    BB = (B * B) % P
    C = (AA - BB) % P
    new_X0 = (AA * BB) % P
    new_Z0 = (C * (AA + (A24 * C) % P)) % P

    # --- Differential addition of R0 and R1 ---
    D = (X1 + Z1) % P
    E = (X1 - Z1) % P
    DA = (E * A) % P
    CB = (D * B) % P
    new_X1 = ((DA + CB) % P) ** 2 % P
    new_Z1 = (xP * (((DA - CB) % P) ** 2)) % P

    return new_X0, new_Z0, new_X1, new_Z1

# --- Encoding / Decoding Helpers ---

def decode_u(u_bytes: bytes) -> int:
    """
    Decode a 32-byte string to an integer,
    masking off the high bit as per RFC 7748.
    """
    u = bytearray(u_bytes)
    u[31] &= 127
    return int.from_bytes(u, "little")

def encode_u(u_int: int) -> bytes:
    """
    Encode an integer as a 32-byte little-endian string.
    """
    return u_int.to_bytes(32, "little")

# --- Montgomery Ladder Implementation ---

def x25519_ladder(k_int: int, u_int: int) -> int:
    """
    Compute the X25519 scalar multiplication using the Montgomery ladder.
    
    This function uses the abstractions for affine/projective conversion
    and a separate ladder step. The base point's affine x-coordinate is u_int.
    
    The two projective points used in the ladder are:
        R0 = (X0:Z0) and R1 = (X1:Z1)
    They are initialized as:
        R0 = (1:0)   (the point-at-infinity in this coordinate system)
        R1 = (u_int:1)  (the affine base point in projective form)
    
    For each bit (from bit 254 down to 0) of the clamped scalar k_int,
    a conditional swap is performed followed by a ladder step.
    
    Finally, the resulting projective coordinate R0 is converted back to affine.
    """
    xP = u_int  # The base point's affine x-coordinate

    # Initialize projective points: R0 = (1:0), R1 = (xP:1)
    X0, Z0 = 1, 0
    X1, Z1 = xP, 1

    swap = 0
    # Process 255 bits (bit 254 down to 0, since the scalar is clamped)
    for t in reversed(range(255)):
        k_t = (k_int >> t) & 1
        swap_bit = swap ^ k_t

        # Conditionally swap the points based on the current bit.
        X0, X1 = cswap(swap_bit, X0, X1)
        Z0, Z1 = cswap(swap_bit, Z0, Z1)
        swap = k_t

        # Perform one simultaneous doubling and differential addition.
        X0, Z0, X1, Z1 = ladder_step(X0, Z0, X1, Z1, xP)

    # Final conditional swap.
    X0, X1 = cswap(swap, X0, X1)
    Z0, Z1 = cswap(swap, Z0, Z1)

    # Convert R0 from projective to affine coordinates.
    return projective_to_affine(X0, Z0)

def x25519(k_bytes: bytes, u_bytes: bytes) -> bytes:
    """
    Public API for X25519 scalar multiplication:
      - Clamps the 32-byte scalar.
      - Decodes the base point from u_bytes.
      - Performs the Montgomery ladder.
      - Returns the result as a 32-byte little-endian string.
    """
    k_int = clamp_scalar(bytearray(k_bytes))
    u_int = decode_u(u_bytes)
    result_int = x25519_ladder(k_int, u_int)
    return encode_u(result_int)


# --- TEST / DEMO CODE BELOW ---
if __name__ == "__main__":
    # Test vector from RFC 7748.
    k_hex = (
        "a546e36bf0527c9d3b16154b82465edd"
        "62144c0ac1fc5a18506a2244ba449ac4"
    )
    u_hex = (
        "e6db6867583030db3594c1a424b15f7c"
        "726624ec26b3353b10a903a6d0ab1c4c"
    )
    expect_hex = (
        "c3da55379de9c6908e94ea4df28d084f"
        "32eccf03491c71f754b4075577a28552"
    )

    k_bytes_ = bytes.fromhex(k_hex)
    u_bytes_ = bytes.fromhex(u_hex)
    result = x25519(k_bytes_, u_bytes_)
    print("Test vector result:", result.hex())
    assert expect_hex == result.hex()

    # Also demonstrate computing a public key from a random private scalar.
    basepoint = b"\x09" + b"\x00" * 31
    import os
    private_scalar = os.urandom(32)  # In a real application use secure randomness.
    public_u = x25519(private_scalar, basepoint)
    print("Random public key:", public_u.hex())
