from abc import ABC, abstractmethod
from util import clamp_scalar, cswap, decode_u, encode_u_coordinate, modinv, projective_to_affine
"""
A simple, from-scratch X25519 implementation in Python using only
built-in big integer support and the Montgomery ladder from RFC 7748.

WARNING: This code is not constant-time and is unsuitable for
production use. It is intended only as a pedagogical reference.
"""

class MontgomeryLadder(ABC):
    def __init__(self):
       self.p = 2**255 - 19
       self.a24 = 121665

    def __call__(self, k_bytes: bytes, u_bytes: bytes) -> bytes:
        return self.x25519(k_bytes, u_bytes)  
       
    @abstractmethod
    def scalar_mult(self, k_int: int, u_int: int) -> int:
        raise NotImplementedError

    def x25519(self, k_bytes: bytes, u_bytes: bytes) -> bytes:
      """
      Perform X25519 scalar multiplication using the optimized Montgomery ladder.

      Steps:
        1. Clamp the 32-byte scalar (per RFC 7748).
        2. Decode the base points x-coordinate from u_bytes.
        3. Run the optimized ladder to compute the scalar multiple.
        4. Encode the resulting x-coordinate as 32 bytes.
      """
      k_int = clamp_scalar(bytearray(k_bytes))
      xP = decode_u(u_bytes)
      result_int = self.scalar_mult(k_int, xP)
      return encode_u_coordinate(result_int)
    

class MontgomeryLadderRFC7748(MontgomeryLadder):

    def scalar_mult(self, k_int: int, u_int: int) -> int:
        """
        Perform the Montgomery ladder (scalar multiplication) on Curve25519:
        X25519(k, u) = k * (u : 1) in the group law, returning the
        resulting u-coordinate as an integer.

        Follows the pseudo-code in RFC 7748, section 5.
        """
        x1 = u_int
        x2, z2 = 1, 0
        x3, z3 = u_int, 1
        swap = 0

        # Loop over bits of k from top (254) down to 0
        # (Note the top bit is always set after clamping, so we skip bit 255)

        # TODO: this works, but we can do better, e.g. look at Martin's tutorial
        # this is basically from the RF
        for t in reversed(range(255)):
            k_t = (k_int >> t) & 1
            if k_t != swap:
                x2, x3 = x3, x2
                z2, z3 = z3, z2
                swap = k_t

            # The curve arithmetic
            #  -- all operations mod P
            A = (x2 + z2) % self.p
            B = (x2 - z2) % self.p
            AA = (A * A) % self.p
            BB = (B * B) % self.p
            E = (AA - BB) % self.p
            C = (x3 + z3) % self.p
            D = (x3 - z3) % self.p
            DA = (D * A) % self.p    
            CB = (C * B) % self.p
            x3 = (DA + CB) % self.p
            x3 = (x3 * x3) % self.p
            z3 = (DA - CB) % self.p
            z3 = (z3 * z3) % self.p
            z3 = (z3 * x1) % self.p
            x2 = (AA * BB) % self.p
            z2 = (E * ((AA + (self.a24 * E) % self.p) % self.p)) % self.p

        # Last swap if needed
        if swap:
            x2, x3 = x3, x2
            z2, z3 = z3, z2

        # Return x2 * (z2^(p-2)) mod p  (the u-coordinate)
        # We use Fermat's little theorem for the inverse: z2^(p-2) mod p
        # inv_z2 = pow(z2, P - 2, P)
        inv_z2: int = modinv(z2, self.p)
        return (x2 * inv_z2) % self.p


class MontgomeryLadderOptimized(MontgomeryLadder):

    def scalar_mult(self, k_int: int, xP: int) -> int:
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
            a, b = cswap(bit, a, b, self.p)
            c, d = cswap(bit, c, d, self.p)

            # --- Begin ladder step ---
            # v1 = a + c
            e = a + c % self.p
            # v2 = a - c
            a = a - c % self.p
            # v3 = b + d
            c = b + d % self.p
            # v4 = b - d
            b = b - d % self.p
            # v5 = (v1)^2 = (a+c)^2
            d = e * e % self.p
            # v6 = (v2)^2 = (a-c)^2
            f_val = a * a % self.p
            # v7 = (b + d)* (a - c)
            a = c * a % self.p
            # v8 = (b - d) * (a + c)
            c = b * e % self.p
            # v9 = v7 + v8
            e = a + c % self.p
            # v10 = v7 - v8
            a = a - c % self.p
            # v11 = (v10)^2
            b = a * a % self.p
            # v12 = v5 - v6
            c = d - f_val % self.p
            # v13 = 121665 * v12
            a = c * self.a24 % self.p
            # v14 = v13 + v5
            a = a + d % self.p
            # v15 = v12 * v14
            c = c * a % self.p
            # v16 = v5 * v6
            a = d * f_val % self.p
            # v17 = v11 * xP
            d = b * xP % self.p
            # v18 = (v9)^2
            b = e * e % self.p
            # --- Final swap (if bit==1) ---
            a, b = cswap(bit, a, b, self.p)
            c, d = cswap(bit, c, d, self.p)
            # End of one ladder iteration

        # After the loop, a and c are the numerator and denominator.
        inv_c = modinv(c, self.p)
        return a * inv_c % self.p
    


class MontgomeryLadderMKTutorial(MontgomeryLadder):
    # This is 10 multiplications
    def ladder_step(self, X0: int, Z0: int, X1: int, Z1: int, xP: int) -> tuple[int, int, int, int]:
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
        A = (X0 + Z0) % self.p
        B = (X0 - Z0) % self.p
        AA = (A * A) % self.p
        BB = (B * B) % self.p
        C = (AA - BB) % self.p
        new_X0 = (AA * BB) % self.p
        new_Z0 = (C * (AA + (self.a24 * C) % self.p)) % self.p

        # --- Differential addition of R0 and R1 ---
        D = (X1 + Z1) % self.p
        E = (X1 - Z1) % self.p
        DA = (E * A) % self.p
        CB = (D * B) % self.p
        new_X1 = ((DA + CB) % self.p) ** 2 % self.p
        new_Z1 = (xP * (((DA - CB) % self.p) ** 2)) % self.p

        return new_X0, new_Z0, new_X1, new_Z1

    def scalar_mult(self, k_int: int, u_int: int) -> int:
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
            X0, X1 = cswap(swap_bit, X0, X1, self.p)
            Z0, Z1 = cswap(swap_bit, Z0, Z1, self.p)
            swap = k_t

            # Perform one simultaneous doubling and differential addition.
            X0, Z0, X1, Z1 = self.ladder_step(X0, Z0, X1, Z1, xP)

        # Final conditional swap.
        X0, X1 = cswap(swap, X0, X1, self.p)
        Z0, Z1 = cswap(swap, Z0, Z1, self.p)

        # Convert R0 from projective to affine coordinates.
        return projective_to_affine(X0, Z0, self.p)