import unittest

from parameterized import parameterized

from curve import IdentityPoint
from ed25519.affine_edwards_curve import AffineEdwardsCurve
from ed25519.extended_edwards_curve import ExtendedEdwardsCurve


class TestEdwardsCurveOperations(unittest.TestCase):
    @parameterized.expand(
        [
            (AffineEdwardsCurve(), "AffineEdwardsCurve"),
            (ExtendedEdwardsCurve(), "ExtendedEdwardsCurve"),
        ]
    )  # type: ignore
    def test_compress_uncompress(self, curve, curve_name) -> None:
        """
        Test that compressing a non-identity point and then uncompressing it
        yields a point equal to the original (via point_equals).
        """
        # Use the curve’s base point.
        P = curve.B
        compressed = curve.compress(P)
        self.assertEqual(
            len(compressed),
            32,
            f"Compressed point length is not 32 for {curve_name}",
        )

        uncompressed = curve.uncompress(compressed)
        self.assertTrue(
            curve.point_equals(P, uncompressed),
            f"Uncompressed point does not equal original for {curve_name}",
        )

    @parameterized.expand(
        [
            (AffineEdwardsCurve(), "AffineEdwardsCurve"),
            (ExtendedEdwardsCurve(), "ExtendedEdwardsCurve"),
        ]
    )  # type: ignore
    def test_compress_identity_raises(self, curve, curve_name) -> None:
        """Test that attempting to compress the identity element raises a ValueError."""
        with self.assertRaises(
            ValueError, msg=f"Compressing identity did not raise error for {curve_name}"
        ):
            _ = curve.compress(IdentityPoint)

    @parameterized.expand(
        [
            (AffineEdwardsCurve(), "AffineEdwardsCurve"),
            (ExtendedEdwardsCurve(), "ExtendedEdwardsCurve"),
        ]
    )  # type: ignore
    def test_uncompress_invalid_length(self, curve, curve_name) -> None:
        """
        Test that providing a byte string of length not equal to 32 to uncompress
        raises a ValueError.
        """
        invalid_bytes = b"\x00" * 31  # Only 31 bytes instead of 32.
        with self.assertRaises(
            ValueError,
            msg=f"Uncompressing invalid length did not raise error for {curve_name}",
        ):
            _ = curve.uncompress(invalid_bytes)

    @parameterized.expand(
        [
            (AffineEdwardsCurve(), "AffineEdwardsCurve"),
            (ExtendedEdwardsCurve(), "ExtendedEdwardsCurve"),
        ]
    )  # type: ignore
    def test_point_equals_identity(self, curve, curve_name) -> None:
        """
        Test that the point_equals method treats IdentityPoint correctly.
        (Note: per our implementation, if either point is IdentityPoint,
         point_equals returns True.)
        """
        P = curve.B
        self.assertTrue(
            curve.point_equals(P, IdentityPoint),
            f"Point equality with IdentityPoint failed for {curve_name}",
        )
        self.assertTrue(
            curve.point_equals(IdentityPoint, P),
            f"Point equality with IdentityPoint failed for {curve_name}",
        )
        self.assertTrue(
            curve.point_equals(IdentityPoint, IdentityPoint),
            f"IdentityPoint not equal to itself for {curve_name}",
        )

    @parameterized.expand(
        [
            (AffineEdwardsCurve(), "AffineEdwardsCurve"),
            (ExtendedEdwardsCurve(), "ExtendedEdwardsCurve"),
        ]
    )  # type: ignore
    def test_add_identity(self, curve, curve_name) -> None:
        """Test that adding IdentityPoint returns the other point unchanged."""
        P = curve.B
        sum1 = curve.add(P, IdentityPoint)
        sum2 = curve.add(IdentityPoint, P)
        self.assertTrue(
            curve.point_equals(sum1, P),
            f"Adding IdentityPoint did not yield the (P + I) for {curve_name}",
        )
        self.assertTrue(
            curve.point_equals(sum2, P),
            f"Adding IdentityPoint did not yield the (I + P) for {curve_name}",
        )

    @parameterized.expand(
        [
            (AffineEdwardsCurve(), "AffineEdwardsCurve"),
            (ExtendedEdwardsCurve(), "ExtendedEdwardsCurve"),
        ]
    )  # type: ignore
    def test_double_identity(self, curve, curve_name) -> None:
        """Test that doubling IdentityPoint returns IdentityPoint."""
        doubled = curve.double(IdentityPoint)
        self.assertEqual(
            doubled,
            IdentityPoint,
            f"Doubling IdentityPoint did not yield IdentityPoint for {curve_name}",
        )

    @parameterized.expand(
        [
            (AffineEdwardsCurve(), "AffineEdwardsCurve"),
            (ExtendedEdwardsCurve(), "ExtendedEdwardsCurve"),
        ]
    )  # type: ignore
    def test_add_commutativity(self, curve, curve_name) -> None:
        """
        Test that point addition is commutative.
        (P + Q should equal Q + P.)
        """
        P = curve.B
        # Obtain a second point by doubling the base point.
        Q = curve.double(P)
        sum1 = curve.add(P, Q)
        sum2 = curve.add(Q, P)
        self.assertTrue(
            curve.point_equals(sum1, sum2),
            f"Point addition is not commutative for {curve_name}",
        )

    @parameterized.expand(
        [
            (AffineEdwardsCurve(), "AffineEdwardsCurve"),
            (ExtendedEdwardsCurve(), "ExtendedEdwardsCurve"),
        ]
    )  # type: ignore
    def test_double_consistency(self, curve, curve_name) -> None:
        """
        Test that doubling a point is consistent with adding the point to itself.
        (i.e. double(P) should equal add(P, P).)
        """
        P = curve.B
        doubled = curve.double(P)
        added = curve.add(P, P)
        self.assertTrue(
            curve.point_equals(doubled, added),
            f"Doubling and adding a point to itself are not consistent for {curve_name}",
        )


if __name__ == "__main__":
    unittest.main()
