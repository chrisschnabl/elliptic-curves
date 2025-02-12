    
from typing import override
from curve import AffinePoint, Point, IdentityPoint, DoubleAndAddCurve
from tonelli_shanks import tonelli
from util import modinv

from x25519.x25519_curve import X25519Curve


class X25519CurveGroupLaw(X25519Curve, DoubleAndAddCurve):
    def __init__(self):
        self.p = 2**255 - 19
        self.A = 486662
        # Standard base point for Curve25519 (full coordinates).
        self.B: Point = AffinePoint(
            9,
            #14781619447589544791020593568409986887264606134616475288964881837755586237401,
            43114425171068552920764898935933967039370386198203806730763910166200978582548,
        )

    @override
    def decode_public_key(self, public_key_bytes: bytes) -> int:
        if len(public_key_bytes) != 32:
            raise ValueError("Public key must be 32 bytes long")
        # Create a mutable copy and clear bit 255 of the last byte.
        pk_bytes = bytearray(public_key_bytes)
        pk_bytes[31] &= 0x7F  # mask with 0x7F to clear the top bit
        # Convert the little-endian byte array to an integer.
        return int.from_bytes(pk_bytes, byteorder="little")
    
    @override
    def recover_point(self, x: int) -> Point:
        """
        Given an x-coordinate, recover a point (x,y) on the curve defined by:
                y² = x³ + A*x² + x   (mod p)
        by computing a square root of the right-hand side.
        (If there is no square root, a ValueError is raised.)
        Note: In our variant the full point (x and y) is transmitted, so this function is not used.
        """
        rhs = (pow(x, 3, self.p) + self.A * pow(x, 2, self.p) + x) % self.p
        y = tonelli(rhs, self.p)
        if y is None:
            raise ValueError("No valid y for given x")
        
        # Choose the smaller square root, but should be the same
        if y > self.p - y:
            y = self.p - y
        return AffinePoint(x, y)


    @override
    def add(self, P: Point, Q: Point) -> Point:
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
        if P is IdentityPoint:
            return Q
        if Q is IdentityPoint:
            return P

        x1, y1 = P.x, P.y
        x2, y2 = Q.x, Q.y

        if x1 == x2:
            # If y1 + y2 = 0 mod p, then Q is the inverse of P.
            if (y1 + y2) % self.p == 0:
                return IdentityPoint
            # Otherwise, P == Q and we perform doubling.
            return self.double(P)

        lam = ((y2 - y1) * modinv(x2 - x1, self.p)) % self.p
        x3 = (lam * lam - self.A - x1 - x2) % self.p
        y3 = (lam * (x1 - x3) - y1) % self.p
        return AffinePoint(x3, y3)

    @override
    def double(self, P: Point) -> Point:
        """
        Double the point P on the Montgomery curve.
        Uses the tangent line at P:
            λ = (3*x² + 2*A*x + 1) / (2*y) mod p,
            x₃ = λ² - A - 2*x mod p,
            y₃ = λ*(x - x₃) - y mod p.
        Returns the doubled point, or None if y == 0.
        """
        if P is IdentityPoint:
            return IdentityPoint
        x, y = P.x, P.y
        if y == 0:
            return IdentityPoint
        lam = ((3 * x * x + 2 * self.A * x + 1) * modinv(2 * y, self.p)) % self.p
        x3 = (lam * lam - self.A - 2 * x) % self.p
        y3 = (lam * (x - x3) - y) % self.p
        return AffinePoint(x3, y3)

                