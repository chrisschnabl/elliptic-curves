from binascii import unhexlify
import unittest

from util import encode_u_coordinate
from x25519.double_and import BasePoint, point_add, point_double, recover_point, scalar_mult

class TestMontgomeryCurveOperations(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        self.p = 2**255 - 19
        self.A = 486662
        super().__init__(*args, **kwargs)

    def test_recover_point_base(self):
        """
        Test that given the base x-coordinate (9), recover_point returns the
        full affine point (x,y) as specified by RFC 7748.
        """
        P = recover_point(9)
        self.assertEqual(
            P, 
            BasePoint, 
            "recover_point(9) did not return the expected BasePoint."
        )
    
    def test_point_lies_on_curve(self):
        """
        Verify that a recovered point (or one produced by our group operations)
        satisfies the curve equation: y^2 = x^3 + A*x^2 + x  mod p.
        """
        # Test on the base point.
        x, y = recover_point(9)
        lhs = (y * y) % self.p
        rhs = (pow(x, 3, self.p) + self.A * pow(x, 2, self.p) + x) % self.p
        self.assertEqual(lhs, rhs, "The recovered point does not lie on the curve.")
    
    def test_point_double(self):
        """
        Test that doubling the BasePoint yields a point that lies on the curve.
        (Also indirectly tests the modular inversion and tonelli routines.)
        """
        D = point_double(BasePoint)
        # We expect a non-infinity result for the base point doubling.
        self.assertIsNotNone(D, "Doubling BasePoint returned the point at infinity unexpectedly.")
        x, y = D
        lhs = (y * y) % self.p
        rhs = (pow(x, 3, self.p) + self.A * pow(x, 2, self.p) + x) % self.p
        self.assertEqual(lhs, rhs, "The doubled point does not lie on the curve.")
    
    def test_point_add(self):
        """
        Test that point addition works correctly.
        Since adding the BasePoint to itself should be equivalent to doubling it,
        we verify that point_add(BasePoint, BasePoint) equals point_double(BasePoint).
        """
        S = point_add(BasePoint, BasePoint)
        D = point_double(BasePoint)
        self.assertEqual(S, D, "Point addition (P + P) does not equal point doubling.")
    
    def test_scalar_mult_basic(self):
        """
        Test scalar multiplication for small scalars:
         - 1 * P should equal P.
         - 2 * P should equal point_double(P).
         Also check that a multiple (e.g. 3*P) yields a point on the curve.
        """
        P = BasePoint
        self.assertEqual(scalar_mult(1, P), P, "1*P does not equal P.")
        self.assertEqual(scalar_mult(2, P), point_double(P), "2*P does not equal point_double(P).")
        
        P3 = scalar_mult(3, P)
        x, y = P3
        lhs = (y * y) % self.p
        rhs = (pow(x, 3, self.p) + self.A * pow(x, 2, self.p) + x) % self.p
        self.assertEqual(lhs, rhs, "3*P does not lie on the curve.")