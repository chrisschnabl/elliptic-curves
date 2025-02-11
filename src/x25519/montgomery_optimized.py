
# -------------------------------------------------------------------

from util import cswap, modinv
from x25519.ecdh import MontgomeryLadder


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