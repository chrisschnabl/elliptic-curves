from typing import override

from curve import AffinePoint, DoubleAndAddCurve, IdentityPoint, Point
from tonelli import tonelli
from util import modinv
from x25519.x25519_curve import X25519Curve


class X25519CurveGroupLaw(X25519Curve, DoubleAndAddCurve):  # type: ignore
    def __init__(self) -> None:
        super().__init__()
        self.A = 486662

    @override  # type: ignore
    def decode_public_key_bytes(self, public_key_bytes: bytes) -> int:
        # Create a mutable copy and clear bit 255 of the last byte.
        pk_bytes = bytearray(public_key_bytes)
        pk_bytes[31] &= 0x7F  # mask with 0x7F to clear the top bit
        # Convert the little-endian byte array to an integer.
        return int.from_bytes(pk_bytes, byteorder="little")

    @override  # type: ignore
    def recover_point(self, x: int) -> Point:
        """
        Given an x-coordinate, recover a point (x,y) on the curve defined by:

        y**2 = x**3 + A*x**2 + x   (mod p) by computing a square root of the RHS.
        If there is no square root, raise a ValueError.
        """
        rhs = (pow(x, 3, self.p) + self.A * pow(x, 2, self.p) + x) % self.p
        y = tonelli(rhs, self.p)
        if y is None:
            raise ValueError("No valid y for given x")

        # Choose the smaller square root, but should be the same
        if y < self.p - y:
            y = self.p - y
        return AffinePoint(x, y)

    @override  # type: ignore
    def add(self, P: Point, Q: Point) -> Point:
        """
        Add two points P and Q (affine coordinates) on the
        Montgomery curve: y**2 = x**3 + A*x*x + x   (mod p)

        For distinct points:
            λ = (y_2 - y_1) / (x_2 - x_1) mod p,
            x₃ = λ*λ - A - x_1 - x_2 mod p,
            y₃ = λ*(x_1 - x_3) - y_1 mod p.
        If P == Q then doubling is performed.
        Returns the resulting point, or None if the result is the point at infinity.
        """
        if P is IdentityPoint or Q is IdentityPoint:
            return Q if P is IdentityPoint else P

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

    @override  # type: ignore
    def double(self, P: Point) -> Point:
        """
        Double the point P on the Montgomery curve.
        Uses the tangent line at P:
            λ = (3*x*x + 2*A*x + 1) / (2*y) mod p,
            x_3 = λ*λ - A - 2*x mod p,
            y_3 = λ*(x - x_3) - y mod p.
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
