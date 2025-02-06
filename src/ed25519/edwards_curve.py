from util import modinv, sqrt_mod

Point = tuple[int, int]
EdwardsPoint = Point | None
Identity = None


class EdwardsCurve:
    """
    An abstraction of a twisted Edwards curve.

    The curve is defined by:
         -x^2 + y^2 = 1 + d*x^2*y^2  (mod p)
    where d is a constant, and the field is F_p.

    Attributes:
      - d: the curve constant.
      - p: the field prime.
      - q: the subgroup order.
      - B: the base point (in affine coordinates as a tuple (x, y)).
    """

    def __init__(self) -> None:
        # Field prime for Ed25519.
        self.p = 2**255 - 19

        # The Edwards curve constant for Ed25519.
        self.d = (-121665 * modinv(121666, self.p)) % self.p

        # Subgroup order for Ed25519.
        self.q = 2**252 + 27742317777372353535851937790883648493

        # Base point B (as specified in RFC8032) in affine coordinates.
        self.B = (
            15112221349535400772501151409588531511454012693041857206046113283949847762202,
            46316835694926478169428394003475163141307993866256225615783033603165251855960,
        )

    def add(self, P: EdwardsPoint, Q: EdwardsPoint) -> EdwardsPoint:
        """
        Add two points P = (x1, y1) and Q = (x2, y2) on the Edwards curve.

        The formulas used are:
            x3 = (x1*y2 + y1*x2) / (1 + d*x1*x2*y1*y2)
            y3 = (y1*y2 + x1*x2) / (1 - d*x1*x2*y1*y2)
        (Division is implemented via multiplication by the modular inverse.)
        """
        if P is Identity:
            return Q
        if Q is Identity:
            return P

        x1, y1 = P
        x2, y2 = Q
        p = self.p
        d = self.d

        denom = d * x1 * x2 * y1 * y2
        inv_denom_x = modinv((1 + denom) % p, p)
        inv_denom_y = modinv((1 - denom) % p, p)
        x3 = ((x1 * y2 + x2 * y1) * inv_denom_x) % p
        y3 = ((x1 * x2 + y1 * y2) * inv_denom_y) % p
        return (x3, y3)

    def scalar_mult(self, scalar: int, P: EdwardsPoint = Identity) -> EdwardsPoint:
        """
        Compute the scalar multiplication [scalar]*P using the double-and-add method.

        (This implementation is not constant time.)
        """
        # TODO CS: do this scalar product in projective coordinates
        # (extended homogenous, RFC 8032
        Q = Identity  # Here, None is the identity element.
        R = P
        while scalar:
            if scalar & 1:
                Q = R if Q is Identity else self.add(Q, R)
            R = self.add(R, R)
            scalar //= 2
        return Q

    def compress(self, P: EdwardsPoint) -> bytes:
        """
        Compress a point P = (x, y) into a 32-byte string using the Ed25519 convention:

        - Encode y as 32 little endian bytes.
        - Set the most-significant bit (of the last byte) to the lsb of x.
        """
        if P is Identity:
            raise ValueError("Cannot compress Identity Element")
        x, y = P
        y_bytes = y.to_bytes(32, "little")
        y_arr = bytearray(y_bytes)
        if x & 1:  # if odd, y is positive
            y_arr[31] |= 0x80
        else:  # if even, y is negative
            y_arr[31] &= 0x7F
        return bytes(y_arr)

    def uncompress(self, comp: bytes) -> EdwardsPoint:
        # TODO CS: look at the slide again
        """
        Uncompress a 32-byte point into its (x, y) coordinates.

        The encoding assumes:
         - The lower 255 bits encode y (little–endian).
         - The most–significant bit holds the least–significant bit of x.

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
        comp_arr[31] &= 0x7F  # Remove the sign bit

        dx = int.from_bytes(comp_arr, "little")
        if dx >= self.p:
            raise ValueError("Decoded y is not in field range")
        dx_squared = (dx * dx) % self.p
        u = (dx_squared - 1) % self.p  # TODO CS: shouldnt this be 1 - dx
        v = (self.d * dx_squared + 1) % self.p  # v = 1 + d*y^2 mod p
        x_sq = (u * modinv(v, self.p)) % self.p
        x = sqrt_mod(x_sq, self.p)
        if (x & 1) != sign:
            x = (-x) % self.p
        return (x, dx)
