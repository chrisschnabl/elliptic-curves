
from typing import override
from ed255192.edwards_curve import ExtendedPoint, Point, IdentityPoint
from ed255192.affine_edwards_curve import AffineEdwardsCurve, AffinePoint
# TODO: am unhappy with the naming

class ExtendedEdwardsCurve(AffineEdwardsCurve):
    """
    Point represented in extended homogeneous coordinates.
    """        
    def __init__(self) -> None:
        super().__init__()
        # Base point B (as specified in RFC8032) in affine coordinates.
        x = 15112221349535400772501151409588531511454012693041857206046113283949847762202
        y = 46316835694926478169428394003475163141307993866256225615783033603165251855960
        self.B = self._from_affine(AffinePoint(x, y))

    def add(
        self,
        P: Point,
        Q: Point
    ) -> Point:
        """
        Add two points P and Q on the curve using extended homogeneous coordinates.
        """
        if P is IdentityPoint or Q is IdentityPoint:
            return P if P is IdentityPoint else Q

        X1, Y1, Z1, T1 = P.x, P.y, P.z, P.t
        X2, Y2, Z2, T2 = Q.x, Q.y, Q.z, Q.t

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

    @override
    def double(self, P: Point) -> Point:
        """
        Double a point P on the curve using extended homogeneous coordinates.
        """
        if P is IdentityPoint:
            return IdentityPoint

        X1, Y1, Z1 = P.x, P.y, P.z
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

    def _from_affine(self, point: Point) -> Point:  # TOOD CS: not quite happy with typing here
        """
        Convert a point from affine coordinates to extended homogeneous coordinates.
        """
        if point is IdentityPoint or type(point) is ExtendedPoint:
            return point

        return ExtendedPoint(point.x, point.y, 1, point.x * point.y % self.p)

    def _to_affine(self, P: ExtendedPoint) -> AffinePoint:
        """
        Convert a point from extended homogeneous coordinates to affine coordinates.
        """
        if P is IdentityPoint:
            return IdentityPoint

        X, Y, Z, _ = P.x, P.y, P.z, P.t
        x = X * pow(Z, -1, self.p) % self.p
        y = Y * pow(Z, -1, self.p) % self.p
        return AffinePoint(x, y)
    
    @override
    def compress(self, point: Point) -> bytes:
        return super().compress(self._to_affine(point))
    
    @override
    def uncompress(self, comp: bytes) -> "ExtendedPoint":
        return self._from_affine(super().uncompress(comp))
    
    @override
    def point_equals(self, P: Point, Q: Point) -> bool:
        return super().point_equals(self._to_affine(P), self._to_affine(Q))