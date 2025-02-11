from util import modinv
from x25519.ecdh import MontgomeryLadder

"""
A simple, from-scratch X25519 implementation in Python using only
built-in big integer support and the Montgomery ladder from RFC 7748.

WARNING: This code is not constant-time and is unsuitable for
production use. It is intended only as a pedagogical reference.
"""


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