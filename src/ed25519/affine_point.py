from typing import Optional
from ed25519.edwards_curve import EdwardsPoint
from util import modinv, sqrt_mod

class AffinePoint(Ed):
    """
    Point represented in affine coordinates.
    Carries a reference to the underlying curve.
    Arithmetic is implemented via conversion to extended coordinates.
    """
    def __init__(self, x: int, y: int) -> None:
        self.x = x
        self.y = y

    def add(self, other: "AffinePoint") -> "AffinePoint":
        """
        Add two points P = (x1, y1) and Q = (x2, y2) on the Edwards curve.
        The formulas used are:
            x3 = (x1*y2 + y1*x2) / (1 + d*x1*x2*y1*y2)
            y3 = (y1*y2 + x1*x2) / (1 - d*x1*x2*y1*y2)
        """
        if self.is_identity(): 
            return other
        if other.is_identity():
            return self

        x1, y1 = self.x, self.y
        x2, y2 = other.x, other.y
        p = self.p
        d = self.d

        denom = d * x1 * x2 * y1 * y2
        inv_denom_x = modinv((1 + denom) % p, p)
        inv_denom_y = modinv((1 - denom) % p, p)
        x3 = ((x1 * y2 + x2 * y1) * inv_denom_x) % p
        y3 = ((x1 * x2 + y1 * y2) * inv_denom_y) % p
        return AffinePoint(x3, y3)

    def scalar_mult(self, scalar: int) -> Optional["AffinePoint"]:
        return super().scalar_mult(scalar)

    def double(self) -> "AffinePoint":
        # TODO CS: Implement this in faster
        #x3 = (x1*y1+y1*x1)/(1+d*x1*x1*y1*y1)
        #  y3 = (y1*y1-a*x1*x1)/(1-d*x1*x1*y1*y1)
        # https://www.hyperelliptic.org/EFD/g1p/auto-twisted.html
        if self.is_identity():
            return self

        x1, y1 = self.x, self.y
        p = self.p
        d = self.d

        denom = d * x1 * x1 * y1 * y1
        inv_denom_x = modinv((1 + denom) % p, p)
        inv_denom_y = modinv((1 - denom) % p, p)
        x3 = ((x1 * y1 + y1 * x1) * inv_denom_x) % p
        y3 = ((y1 * y1 - self.a * x1 * x1) * inv_denom_y) % p
        return AffinePoint(x3, y3)
    

    def compress(self) -> bytes:
        """
        Compress a point P = (x, y) into a 32-byte string using the Ed25519 convention:

        - Encode y as 32 little endian bytes.
        - Set the most-significant bit (of the last byte) to the lsb of x.
        """
        if self.is_identity():
            raise ValueError("Cannot compress Identity Element")
        x, y = self.x, self.y
        y_bytes = y.to_bytes(32, "little")
        y_arr = bytearray(y_bytes)
        if x & 1:  # if odd, y is positive
            y_arr[31] |= 0x80
        else:  # if even, y is negative
            y_arr[31] &= 0x7F
        return bytes(y_arr)

    def uncompress(self, comp: bytes) -> "AffinePoint":
        """
        Uncompress a 32-byte point into its (x, y) coordinates.

        The encoding assumes:
         - The lower 255 bits encode y (little-endian).
         - The most-significant bit holds the least-significant bit of x.

        To recover x, use the curve equation:
           -x^2 + y^2 = 1 + d*x^2*y^2.
        Rearranged to:
           x^2 = (y^2 - 1) / (d*y^2 + 1)   (mod p)
        and then take a square root (using sqrt_mod).
        """
        if len(comp) != 32:
            raise ValueError("Compressed point must be 32 bytes")

        comp_arr = bytearray(comp)
        sign = (comp_arr[31] >> 7) & 1  # Recover the sign of y
        comp_arr[31] &= 0x7F  # Reset the sign bit

        dx = int.from_bytes(comp_arr, "little")
        if dx >= self.p:
            raise ValueError("Decoded y is not in field range")
        dx_squared = (dx * dx) % self.p
        u = (dx_squared - 1) % self.p
        v = (self.d * dx_squared + 1) % self.p
        x_sq = (u * modinv(v, self.p)) % self.p
        x = sqrt_mod(x_sq, self.p)
        if (x & 1) != sign:
            x = (-x) % self.p
        return (x, dx)

    def equals(self, other: "AffinePoint") -> bool:
        return self.x == other.x and self.y == other.y

    def __repr__(self) -> str:
        return f"AffinePoint(x={self.x}, y={self.y})"
    