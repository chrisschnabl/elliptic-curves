from util import clamp_scalar, modinv  # Assume these are provided

# -------------------------------------------------------------------
# Global parameters for Curve25519
# -------------------------------------------------------------------
P = 2**255 - 19
# The constant 121665 appears in the optimized formulas.
C_121665 = 121665

# -------------------------------------------------------------------
# Field arithmetic helpers (all operations mod P)
# -------------------------------------------------------------------
def fadd(x: int, y: int) -> int:
    return (x + y) % P

def fsub(x: int, y: int) -> int:
    return (x - y) % P

def fmul(x: int, y: int) -> int:
    return (x * y) % P

# -------------------------------------------------------------------
# Constant-time conditional swap (not truly constant-time in Python)
# -------------------------------------------------------------------
def cswap(swap: int, x: int, y: int) -> (int, int):
    """
    Swap x and y if swap == 1, otherwise leave them unchanged.
    (Uses arithmetic so that both branches perform the same operations.)
    """
    dummy = swap * (x - y)
    x_new = (x - dummy) % P
    y_new = (y + dummy) % P
    return x_new, y_new

# -------------------------------------------------------------------
# Encoding and decoding helpers for 32-byte values
# -------------------------------------------------------------------
def decode_u(u_bytes: bytes) -> int:
    """
    Decode a 32-byte little-endian string into an integer.
    (Mask off the high bit per RFC 7748.)
    """
    u = bytearray(u_bytes)
    u[31] &= 127
    return int.from_bytes(u, "little")

def encode_u(u_int: int) -> bytes:
    """
    Encode an integer as a 32-byte little-endian string.
    """
    return u_int.to_bytes(32, "little")

# -------------------------------------------------------------------
# Optimized Montgomery ladder step (per Listing 5 in the tutorial)
#
# The state is maintained in four variables:
#   a, b, c, d
#
# Their interpretation is as follows (using the notation from the tutorial):
#
#   At the start of the iteration:
#       a = Xᵢ,   b = Xᵢ₊₁,   c = Zᵢ,   d = Zᵢ₊₁
#
# and the base point’s affine x-coordinate is xP.
#
# After one iteration (depending on the bit) the new state is either
#
#       (X₂ᵢ, Z₂ᵢ, 4·X₂ᵢ₊₁, 4·Z₂ᵢ₊₁)
#
# or (with a swap)
#
#       (4·X₂ᵢ₊₁, 4·Z₂ᵢ₊₁, X₂ᵢ₊₂, Z₂ᵢ₊₂)
#
# (The extra factors of 4 cancel later.) Finally, after all iterations the
# x-coordinate is recovered as a/c.
# -------------------------------------------------------------------
def x25519_optimized_ladder(k_int: int, xP: int) -> int:
    """
    Compute scalar multiplication using an optimized Montgomery ladder.

    The initial state is set as:
       a = 1,    b = xP,
       c = 0,    d = 1.
    (This corresponds to representing the “point-at-infinity” as (1:0) and the
    base point in projective coordinates as (xP:1).)
    
    We then loop over the 255 bits (bit 254 down to 0) of the clamped scalar,
    updating the state with 18 arithmetic operations per iteration (10 multiplications).
    
    Finally, we return the affine x-coordinate as a * inv(c) mod P.
    """
    # Initialize state
    a = 1
    b = xP
    c = 0
    d = 1

    # Process bits 254 down to 0
    for i in range(254, -1, -1):
        bit = (k_int >> i) & 1

        # --- Pre-step swap (if bit==1) ---
        a, b = cswap(bit, a, b)
        c, d = cswap(bit, c, d)

        # --- Begin ladder step ---
        # v1 = a + c
        e = fadd(a, c)
        # v2 = a - c
        a = fsub(a, c)
        # v3 = b + d
        c = fadd(b, d)
        # v4 = b - d
        b = fsub(b, d)
        # v5 = (v1)^2 = (a+c)^2
        d = fmul(e, e)
        # v6 = (v2)^2 = (a-c)^2
        f_val = fmul(a, a)
        # v7 = (b + d)* (a - c)
        a = fmul(c, a)
        # v8 = (b - d) * (a + c)
        c = fmul(b, e)
        # v9 = v7 + v8
        e = fadd(a, c)
        # v10 = v7 - v8
        a = fsub(a, c)
        # v11 = (v10)^2
        b = fmul(a, a)
        # v12 = v5 - v6
        c = fsub(d, f_val)
        # v13 = 121665 * v12
        a = fmul(c, C_121665)
        # v14 = v13 + v5
        a = fadd(a, d)
        # v15 = v12 * v14
        c = fmul(c, a)
        # v16 = v5 * v6
        a = fmul(d, f_val)
        # v17 = v11 * xP
        d = fmul(b, xP)
        # v18 = (v9)^2
        b = fmul(e, e)
        # --- Final swap (if bit==1) ---
        a, b = cswap(bit, a, b)
        c, d = cswap(bit, c, d)
        # End of one ladder iteration

    # After the loop, a and c are the numerator and denominator.
    inv_c = modinv(c, P)
    return fmul(a, inv_c)

# -------------------------------------------------------------------
# Public API: X25519 using the optimized ladder.
# -------------------------------------------------------------------
"""
1. Use Projective Coordinates
Projective Representation:
Represent points as 
(
𝑋
:
𝑍
)
(X:Z) so that division (i.e., field inversion) is deferred until after all iterations.
Advantage:
This avoids expensive inversion operations inside the main loop.
2. Merge Doubling and Differential Addition
Combined Formulas:
Instead of computing doubling and differential addition separately, the two operations are merged into a single set of operations.
Benefit:
This allows the reuse of common subexpressions, reducing the overall number of finite field multiplications.
3. Reuse Intermediate Values
Define the following temporary variables:

𝑣
1
=
𝑎
+
𝑐
v 
1
​
 =a+c
𝑣
2
=
𝑎
−
𝑐
v 
2
​
 =a−c
𝑣
3
=
𝑏
+
𝑑
v 
3
​
 =b+d
𝑣
4
=
𝑏
−
𝑑
v 
4
​
 =b−d
Then compute:

𝑣
5
=
𝑣
1
2
v 
5
​
 =v 
1
2
​
 
𝑣
6
=
𝑣
2
2
v 
6
​
 =v 
2
2
​
 
𝑣
7
=
𝑣
3
⋅
𝑣
2
v 
7
​
 =v 
3
​
 ⋅v 
2
​
 
𝑣
8
=
𝑣
4
⋅
𝑣
1
v 
8
​
 =v 
4
​
 ⋅v 
1
​
 
𝑣
9
=
𝑣
7
+
𝑣
8
v 
9
​
 =v 
7
​
 +v 
8
​
 
𝑣
10
=
𝑣
7
−
𝑣
8
v 
10
​
 =v 
7
​
 −v 
8
​
 
These reuse steps ensure that common calculations are performed only once.

4. Reduce Multiplications
Reduction:
With careful combination and reuse of intermediate values, each ladder iteration requires only 10 multiplications (along with additional additions/subtractions) instead of 14.
Significance:
Multiplications are the most costly field operations, so reducing their number has a direct impact on performance.
5. Cancel Extraneous Factors
Extra Factors:
The differential addition branch may introduce extra factors (e.g., factors of 4).
Cancellation:
These factors cancel out when converting the result from projective coordinates back to affine coordinates, so they do not affect the final result.
6. Conditional Swaps for Constant-Time Behavior
Branchless Swaps:
Use constant-time conditional swaps (e.g., a branchless cswap function) to swap the roles of the two projective points based on the current bit of the scalar.
Purpose:
This ensures that the algorithm executes in constant time, preventing timing attacks.

"""



def x25519(k_bytes: bytes, u_bytes: bytes) -> bytes:
    """
    Perform X25519 scalar multiplication using the optimized Montgomery ladder.

    Steps:
      1. Clamp the 32-byte scalar (per RFC 7748).
      2. Decode the base point’s x-coordinate from u_bytes.
      3. Run the optimized ladder to compute the scalar multiple.
      4. Encode the resulting x-coordinate as 32 bytes.
    """
    k_int = clamp_scalar(bytearray(k_bytes))
    xP = decode_u(u_bytes)
    result_int = x25519_optimized_ladder(k_int, xP)
    return encode_u(result_int)

# -------------------------------------------------------------------
# Test / Demo Code
# -------------------------------------------------------------------
if __name__ == "__main__":
    # Test vector from RFC 7748 (should match the expected output)
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
    result = x25519_optimized(k_bytes_, u_bytes_)
    print("Optimized test vector result:", result.hex())
    assert result.hex() == expect_hex

    # Also demonstrate computing a public key from a random scalar.
    basepoint = b"\x09" + b"\x00" * 31
    import os
    private_scalar = os.urandom(32)
    public_key = x25519_optimized(private_scalar, basepoint)
    print("Optimized random public key:", public_key.hex())
