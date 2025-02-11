# -------------------------------------------------------------------
# Constant-time conditional swap (not truly constant-time in Python)
# -------------------------------------------------------------------
def cswap(swap: int, x: int, y: int, p: int) -> tuple[int, int]:
    """
    Swap x and y if swap == 1, otherwise leave them unchanged.
    (Uses arithmetic so that both branches perform the same operations.)
    """
    dummy = swap * (x - y)
    x_new = (x - dummy) % p
    y_new = (y + dummy) % p
    return x_new, y_new

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

def decode_u(u_bytes: bytes) -> int:
    """
    Decode a 32-byte little-endian string into an integer.
    (Mask off the high bit per RFC 7748.)
    """
    u = bytearray(u_bytes)
    u[31] &= 127
    return int.from_bytes(u, "little")

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