import unittest

from x25519.group_law import X25519CurveGroupLaw


class TestMontgomeryCurveOperations(unittest.TestCase):
    def __init__(self, *args, **kwargs) -> None:  # type: ignore
        self.curve = X25519CurveGroupLaw()
        super().__init__(*args, **kwargs)

    def test_recover_point_base(self) -> None:
        """
        Test that given the base x-coordinate (9), recover_point returns the full affine
        point (x,y) as specified by RFC 7748.
        """
        P = self.curve.recover_point(9)
        self.assertEqual(
            P, self.curve.B, "recover_point(9) did not return the expected BasePoint."
        )

    def test_point_lies_on_curve(self) -> None:
        """
        Verify that a recovered point (or one produced by our group operations)
        satisfies the curve equation: y^2 = x^3 + A*x^2 + x  mod p.
        """
        # Test on the base point.
        P = self.curve.recover_point(9)
        x, y = P.x, P.y
        lhs = (y * y) % self.curve.p
        rhs = (
            pow(x, 3, self.curve.p) + self.curve.A * pow(x, 2, self.curve.p) + x
        ) % self.curve.p
        self.assertEqual(lhs, rhs, "The recovered point does not lie on the curve.")

    def test_point_double(self) -> None:
        """
        Test that doubling the BasePoint yields a point that lies on the curve.

        (Also indirectly tests the modular inversion and tonelli routines.)
        """
        D = self.curve.double(self.curve.B)
        # We expect a non-infinity result for the base point doubling.
        self.assertIsNotNone(
            D, "Doubling BasePoint returned the point at infinity unexpectedly."
        )
        x, y = D.x, D.y
        lhs = (y * y) % self.curve.p
        rhs = (
            pow(x, 3, self.curve.p) + self.curve.A * pow(x, 2, self.curve.p) + x
        ) % self.curve.p
        self.assertEqual(lhs, rhs, "The doubled point does not lie on the curve.")

    def test_point_add(self) -> None:
        """
        Test that point addition works correctly.

        Since adding the BasePoint to itself should be equivalent to doubling it, we
        verify that point_add(BasePoint, BasePoint) equals point_double(BasePoint).
        """
        S = self.curve.add(self.curve.B, self.curve.B)
        D = self.curve.double(self.curve.B)
        self.assertEqual(S, D, "Point addition (P + P) does not equal point doubling.")

    def test_scalar_mult_basic(self) -> None:
        """
        Test scalar multiplication for small scalars:
         - 1 * P should equal P.
         - 2 * P should equal point_double(P).
         Also check that a multiple (e.g. 3*P) yields a point on the curve.
        """
        P = self.curve.B
        self.assertEqual(self.curve.scalar_mult(P, 1), P, "1*P does not equal P.")
        self.assertEqual(
            self.curve.scalar_mult(P, 2),
            self.curve.double(P),
            "2*P does not equal point_double(P).",
        )

        P3 = self.curve.scalar_mult(P, 3)
        x, y = P3.x, P3.y
        lhs = (y * y) % self.curve.p
        rhs = (
            pow(x, 3, self.curve.p) + self.curve.A * pow(x, 2, self.curve.p) + x
        ) % self.curve.p
        self.assertEqual(lhs, rhs, "3*P does not lie on the curve.")
