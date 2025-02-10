from typing import Tuple, Optional

from util import modinv, sqrt_mod


class TwistedEdwardsCurve:
    def __init__(self) -> None:
        self.a = -1
        self.p = 2**255 - 19
        # The Edwards curve constant for Ed25519.
        self.d = (-121665 * modinv(121666, self.p)) % self.p
        # The subgroup order for Ed25519.
        self.q = 2**252 + 27742317777372353535851937790883648493
        # Base point B (as specified in RFC8032) in affine coordinates.
        x = 15112221349535400772501151409588531511454012693041857206046113283949847762202
        y = 46316835694926478169428394003475163141307993866256225615783033603165251855960
        self.B = self._from_affine(x, y)

    def add(
        self,
        P: Optional[Tuple[int, int, int, int]],
        Q: Optional[Tuple[int, int, int, int]]
    ) -> Optional[Tuple[int, int, int, int]]:
        """
        Add two points P and Q on the curve using extended homogeneous coordinates.
        """
        if P is None or Q is None:
            return None

        X1, Y1, Z1, T1 = P
        X2, Y2, Z2, T2 = Q

        A = (Y1-X1)*(Y2-X2) % self.p
        B = (Y1+X1)*(Y2+X2) % self.p
        C = 2 * self.d * T1 * T2 % self.p
        D = 2 * Z1 * Z2 % self.p
        
        E = B-A
        F = D-C
        G = D+C
        H = B+A
        X3 = E*F % self.p
        Y3 = G*H % self.p
        T3 = E*H % self.p
        Z3 = F*G % self.p

        return (X3, Y3, Z3, T3)
    
    def double(self, P: Optional[Tuple[int, int, int, int]]) -> Optional[Tuple[int, int, int, int]]:
        """
        Double a point P on the curve using extended homogeneous coordinates.
        """
        #return self.add(P, P)
        if P is None:
            return None

        X1, Y1, Z1, _ = P
        # Formulae 7 from the point doubling paper
        # TODO: CS why don't we use T1?
        A = X1**2 % self.p
        B = Y1**2 % self.p
        C = 2*Z1**2 % self.p
        D = self.a*A % self.p
        E = (X1+Y1)**2 - A - B % self.p
        G = D + B % self.p
        F = G - C % self.p
        H = D - B % self.p
        X3 = E*F % self.p
        Y3 = G*H % self.p
        T3 = E*H % self.p
        Z3 = F*G % self.p

        return (X3, Y3, Z3, T3)

    def scalar_mult(
        self, k: int, P: Tuple[int, int, int, int]
    ) -> Tuple[int, int, int, int]:
        """
        Perform scalar multiplication k * P using the double-and-add method.
        """
        Q = None
        R = P
        while k > 0:
            if k % 2 == 1:
                Q = R if Q is None else self.add(Q, R)
            R = self.double(R)
            k //= 2
        return Q
    
    def _from_affine(self, x: int, y: int) -> Tuple[int, int, int, int]:
        """
        Convert a point from affine coordinates to extended homogeneous coordinates.
        """
        return (x, y, 1, x * y % self.p)

    def _to_affine(self, P: Optional[Tuple[int, int, int, int]]) -> Optional[Tuple[int, int]]:
        """
        Convert a point from extended homogeneous coordinates to affine coordinates.
        """
        if P is None:
            return None

        X, Y, Z, T = P
        x = X * pow(Z, -1, self.p) % self.p
        y = Y * pow(Z, -1, self.p) % self.p
        return (x, y)
    

    def _compress(self, P: Tuple[int, int]) -> bytes:
        """
        Compress a point P = (x, y) into a 32-byte string using the Ed25519 convention:

        - Encode y as 32 little endian bytes.
        - Set the most-significant bit (of the last byte) to the lsb of x.
        """
        x, y = P
        y_bytes = y.to_bytes(32, "little")
        y_arr = bytearray(y_bytes)
        if x & 1:  # if odd, y is positive
            y_arr[31] |= 0x80
        else:  # if even, y is negative
            y_arr[31] &= 0x7F
        return bytes(y_arr)
    
    def compress(self, P: Tuple[int, int, int, int]) -> bytes:
        return self._compress(self._to_affine(P))

    def _uncompress(self, comp: bytes) -> Tuple[int, int]:
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

    def uncompress(self, comp: bytes) -> Tuple[int, int, int, int]:
        """
        Uncompress a 32-byte point into its (x, y) coordinates.
        """
        return self._from_affine(*self._uncompress(comp))
    
    def equals(self, P: Tuple[int, int, int, int], Q: Tuple[int, int, int, int]) -> bool:
        return self._to_affine(P) == self._to_affine(Q)
