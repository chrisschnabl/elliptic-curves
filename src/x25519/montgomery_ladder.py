from curve import AffinePoint, Point
from util import cswap, modinv, projective_to_affine
from x25519.x25519_curve import X25519Curve

"""
WARNING: This code is not constant-time and is unsuitable for
production use. It is intended only as a pedagogical reference.
"""


class MontgomeryLadderRFC7748(X25519Curve):  # type: ignore
    def scalar_mult(self, R: Point, scalar: int) -> Point:
        """
        Perform the Montgomery ladder (scalar multiplication) on Curve25519: X25519(k, u)
        = k * (u : 1) in the group law, returning the resulting u-coordinate as an
        integer.

        Follows the pseudo-code in RFC 7748, section 5.
        """
        u_int = R.x

        x1 = u_int
        x2, z2 = 1, 0
        x3, z3 = u_int, 1
        swap = 0

        # Loop over bits of k from top (254) down to 0
        # (Note the top bit is always set after clamping, so we skip bit 255)
        for t in reversed(range(255)):
            k_t = (scalar >> t) & 1
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
        x = (x2 * inv_z2) % self.p
        return AffinePoint(x, 0)  # We only care about the x-coordinate


class MontgomeryLadderMKTutorial(X25519Curve):  # type: ignore
    def ladder_step(
        self, X0: int, Z0: int, X1: int, Z1: int, xP: int
    ) -> tuple[int, int, int, int]:
        """
        Perform one simultaneous step of the Montgomery ladder:
        This is the step defined in Formula (28) of the tutorial.

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

    def scalar_mult(self, R: Point, scalar: int) -> Point:
        """
        Compute the X25519 scalar multiplication using the Montgomery ladder
        described in Martin Kleppmann's tutorial on elliptic curves.

        This function uses the abstractions for affine/projective conversion
        and a separate ladder step. The base point's affine x-coordinate is u_int.

        The two projective points used in the ladder are:
            R0 = (X0:Z0) and R1 = (X1:Z1)
        They are initialized as:
            R0 = (1:0)   (the point-at-infinity in this coordinate system)
            R1 = (u_int:1)  (the affine base point in projective form)

        For each bit (from bit 254 down to 0) of the clamped scalar k_int,
        a conditional swap is performed followed by a ladder step.

        The resulting projective coordinate R0 is converted back to affine.
        """
        xP = R.x

        # Initialize projective points: R0 = (1:0), R1 = (xP:1)
        X0, Z0 = 1, 0
        X1, Z1 = xP, 1

        swap = 0
        # Process 255 bits (bit 254 down to 0)
        for t in reversed(range(255)):
            k_t = (scalar >> t) & 1
            swap_bit = swap ^ k_t
            X0, X1 = cswap(swap_bit, X0, X1, self.p)
            Z0, Z1 = cswap(swap_bit, Z0, Z1, self.p)
            swap = k_t

            X0, Z0, X1, Z1 = self.ladder_step(X0, Z0, X1, Z1, xP)

        X0, X1 = cswap(swap, X0, X1, self.p)
        Z0, Z1 = cswap(swap, Z0, Z1, self.p)

        x = projective_to_affine(X0, Z0, self.p)
        return AffinePoint(x, 0)


class MontgomeryLadderOptimized(X25519Curve):  # type: ignore
    def ladder_step(
        self, a: int, b: int, c: int, d: int, Rx: int
    ) -> tuple[int, int, int, int]:
        """Perform one step of the Montgomery ladder."""
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
        d = b * Rx % self.p
        # v18 = (v9)^2
        b = e * e % self.p

        return a, b, c, d

    def scalar_mult(self, R: Point, scalar: int) -> Point:
        """
        Compute scalar multiplication using an optimized Montgomery ladder described in
        Listing 5 of Martin Kleppmann's tutorial on elliptic curves.
        https://martin.kleppmann.com/papers/curve25519.pdf

        The initial state is set as:
        a = 1,    b = xP,
        c = 0,    d = 1.
        (This corresponds to representing the "point-at-infinity" as (1:0) and the
        base point in projective coordinates as (xP:1).)

        We then loop over the 255 bits (bit 254 down to 0) of the clamped scalar,
        updating the state with 18 arithmetic operations per iteration (10 muls).

        Finally, return the affine x-coordinate as a * inv(c) mod P.
        """
        # Initialize state
        a = 1
        b = R.x
        c = 0
        d = 1

        # Process bits 254 down to 0
        for i in range(254, -1, -1):
            bit = (scalar >> i) & 1

            # --- Pre-step swap (if bit==1) ---
            a, b = cswap(bit, a, b, self.p)
            c, d = cswap(bit, c, d, self.p)

            a, b, c, d = self.ladder_step(a, b, c, d, R.x)

            # --- Final swap (if bit==1) ---
            a, b = cswap(bit, a, b, self.p)
            c, d = cswap(bit, c, d, self.p)

        # After the loop, a and c are the numerator and denominator.
        inv_c = modinv(c, self.p)
        x = a * inv_c % self.p
        return AffinePoint(x, 0)
