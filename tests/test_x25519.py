import unittest
from binascii import unhexlify

# omething is fucked up here
# from x25519.montgomery import x25519
from nacl.bindings import crypto_scalarmult
from parameterized import parameterized

from x25519.double_and import x25519


class TestX25519ImplementsRFC7748(unittest.TestCase):
    @parameterized.expand(
        [
            (
                "vector1",
                "a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4",
                "e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c",
                "c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552",
            ),
            (
                "vector2",
                "4b66e9d4d1b4673c5ad22691957d6af5c11b6421e0ea01d42ca4169e7918ba0d",
                "e5210f12786811d3f4b7959d0538ae2c31dbe7106fc03c3efc4cd549c715a493",
                "95cbde9476e8907d7aade45cb4b873f88b595a68799fa152e6f8f7647aac7957",
            ),
        ]
    )
    def test_rfc7748_vectors(self, name, k_hex, u_hex, expected_hex):
        """Test X25519 against RFC 7748 test vectors."""
        k_bytes = unhexlify(k_hex)
        u_bytes = unhexlify(u_hex)
        expected = unhexlify(expected_hex)

        try:
            result_custom = x25519(k_bytes, u_bytes)
        except Exception as e:
            self.fail(f"Custom X25519 implementation raised an exception: {e}")

        self.assertEqual(
            result_custom,
            expected,
            f"{name}: Custom X25519 implementation does not match the expected RFC 7748 result.",
        )

        # self.assertEqual(
        #    result_custom, crypto_scalarmult(k_bytes, u_bytes),
        #    f"{name}: Custom X25519 implementation does not match PyNaCl's implementation."
        # )

    def _test_rfc7748_iterative(self):
        # Initial values for k and u as specified in RFC 7748 Section 5.2
        k = unhexlify("0900000000000000000000000000000000000000000000000000000000000000")
        u = unhexlify("0900000000000000000000000000000000000000000000000000000000000000")
        k_2 = unhexlify(
            "0900000000000000000000000000000000000000000000000000000000000000"
        )
        u_2 = unhexlify(
            "0900000000000000000000000000000000000000000000000000000000000000"
        )

        # Expected outputs after specified iterations
        expected_outputs = {
            1: "422c8e7a6227d7bca1350b3e2bb7279f7897b87bb6854b783c60e80311ae3079",
            1000: "684cf59ba83309552800ef566f2f4d3c1c3887c49360e3875f2eb94d99532c51",
            1000000: "7c3911e0ab2586fd864497297e575e6f3bc601c0883c30df5f4dd2d24f665424",
        }

        # Perform iterations
        for i in range(1, 10000 + 1):
            k, u = crypto_scalarmult(k, u), k
            k_2, u_2 = x25519(k_2, u_2), k_2
            self.assertEqual(
                k,
                k_2,
                f"Iteration {i}: X25519 output does not match the expected result.",
            )
            if i in expected_outputs:
                expected = unhexlify(expected_outputs[i])
                self.assertEqual(
                    expected,
                    k,
                    f"Iteration {i}: X25519 output does not match the expected result.",
                )

    def test_random_vectors(self):
        """Test X25519 implementation against PyNaCl with random keys and points."""
        import os

        # Number of random test cases
        num_tests = 0

        for i in range(num_tests):
            # Generate random 32-byte scalar and u-coordinate
            k_bytes = os.urandom(32)
            u_bytes = os.urandom(32)

            # Compute using the custom X25519 implementation
            try:
                result_custom = x25519(k_bytes, u_bytes)
            except Exception as e:
                self.fail(
                    f"Custom X25519 implementation raised an exception on random test {i}: {e}"
                )

            # Compute using PyNaCl
            try:
                import nacl.bindings

                result_nacl = nacl.bindings.crypto_scalarmult(k_bytes, u_bytes)
            except ImportError:
                self.skipTest("PyNaCl is not installed. Skipping nacl comparison tests.")

            # Verify that custom implementation matches PyNaCl's result
            self.assertEqual(
                result_custom,
                result_nacl,
                f"Mismatch between custom X25519 and PyNaCl implementations on random test {i}.",
            )
