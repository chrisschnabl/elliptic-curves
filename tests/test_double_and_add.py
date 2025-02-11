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

    def test_rfc7748_iterative(self):
        """
        Implements the iterative test described in RFC 7748 Section 5.2.
        For each iteration, the procedure is:
          - new_point = scalar_mult(k, u)
          - new k is the x-coordinate of new_point
          - new u is recovered from the old k (since we assume the full point is sent)
          - Then update: k ← new_k, u ← new_u.
        At iterations 1, 1000, and 1,000,000 we check that the encoded x-coordinate
        matches the RFC test vectors.
        """
        # Initial values (RFC 7748 specifies both the scalar and u-coordinate start as 9).
        # Here, we represent k as the integer 9, and u as the full base point recovered from x = 9.
        k = 9
        u = recover_point(9)  # returns the point (9, y) corresponding to the base point

        # Expected outputs (encoded 32-byte little-endian values) after the specified iterations.
        expected_outputs = {
            1: "422c8e7a6227d7bca1350b3e2bb7279f7897b87bb6854b783c60e80311ae3079",
            1000: "684cf59ba83309552800ef566f2f4d3c1c3887c49360e3875f2eb94d99532c51",
            1000000: "7c3911e0ab2586fd864497297e575e6f3bc601c0883c30df5f4dd2d24f665424",
        }

        # Iterate for 1,000,000 iterations.
        for i in range(1, 1000000 + 1):
            # Compute the new point via scalar multiplication: new_point = k * u.
            point = scalar_mult(k, u)  # returns a point (x, y)
            new_k = point[0]           # The new scalar is the x-coordinate of the resulting point.
            # According to the iterative procedure, the new u is the previous k,
            # but since we assume that the full point (including y) is transmitted,
            # we recover it using recover_point.
            new_u = recover_point(k)
            # Update the iterative state.
            k, u = new_k, new_u

            # When we reach an iteration for which an expected output is defined, verify it.
            if i in expected_outputs:
                expected = unhexlify(expected_outputs[i])
                actual = encode_u_coordinate(k)
                self.assertEqual(
                    actual,
                    expected,
                    f"Iteration {i}: X25519 output does not match the expected result."
                )

    def TODO_test_x25519_rfc_test_vector(self):
        """
        Using the RFC 7748 test vector, verify that x25519 produces the expected output.
        
        RFC 7748 (Section 5.2) specifies:
          - Input scalar (k):
              a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4
          - Input u-coordinate (u):
              e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c
          - Expected output u-coordinate:
              c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552
              
        In our implementation, the x25519 function recovers the full point (using y),
        performs scalar multiplication, and then returns the x-coordinate.
        """
        k_hex = "a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4"
        u_hex = "e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c"
        expected_hex = "c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552"
        
        k_bytes = bytes.fromhex(k_hex)
        u_bytes = bytes.fromhex(u_hex)
        result = x25519(k_bytes, u_bytes)
        
        self.assertEqual(
            result, 
            bytes.fromhex(expected_hex),
            "x25519 did not produce the expected output for the RFC test vector."
        )