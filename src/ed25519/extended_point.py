
from typing import override
from ed25519.affine_edwards_curve import AffinePoint
from ed25519.edwards_curve import EdwardsPoint


class ExtendedPoint(EdwardsPoint):
    """
    Point represented in extended homogeneous coordinates.
    """        

    def __init__(self, x: int, y: int, z: int, t: int) -> None:
        super().__init__()
        self.x = x
        self.y = y
        self.z = z
        self.t = t

    def add(
        self,
        other: "ExtendedPoint"
    ) -> "ExtendedPoint":
        """
        Add two points P and Q on the curve using extended homogeneous coordinates.
        """
        if self.is_identity() or other.is_identity():
            return self if self.is_identity() else other

        X1, Y1, Z1, T1 = self.x, self.y, self.z, self.t
        X2, Y2, Z2, T2 = other.x, other.y, other.z, other.t

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

        return ExtendedPoint(X3, Y3, Z3, T3)

        
    def double(self) -> "ExtendedPoint":
        """
        Double a point P on the curve using extended homogeneous coordinates.
        """
        #return self.add(P, P)
        if self.is_identity():
            return self

        X1, Y1, Z1 = self.x, self.y, self.z
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

        return ExtendedPoint(X3, Y3, Z3, T3) 

    def _from_affine(self, point: "AffinePoint") -> None:
        """
        Convert a point from affine coordinates to extended homogeneous coordinates.
        """
        self.x = point.x
        self.y = point.y
        self.z = 1
        self.t = point.x * point.y % self.p

    def _to_affine(self) -> "AffinePoint":
        """
        Convert a point from extended homogeneous coordinates to affine coordinates.
        """
        if self.is_identity():
            return None

        X, Y, Z, _ = self.x, self.y, self.z, self.t
        x = X * pow(Z, -1, self.p) % self.p
        y = Y * pow(Z, -1, self.p) % self.p
        return AffinePoint(x, y)
    
    @override
    def compress(self) -> bytes:
        return self._to_affine().compress()
    
    @override
    def uncompress(self, comp: bytes) -> "ExtendedPoint":
        return self._from_affine(self.uncompress(comp))
    
    def equals(self, other: "ExtendedPoint") -> bool:
        return self._to_affine().equals(other._to_affine())
    
    def __repr__(self) -> str:
        return f"ExtendedPoint(x={self.x}, y={self.y}, z={self.z}, t={self.t})"