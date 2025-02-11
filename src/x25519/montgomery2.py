from util import clamp_scalar, cswap, modinv, projective_to_affine
from x25519.ecdh import MontgomeryLadder


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