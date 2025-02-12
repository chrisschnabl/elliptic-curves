from binascii import unhexlify
from x25519.group_law import X25519CurveGroupLaw
from nacl.bindings import crypto_scalarmult
import unittest


class TestX25519ImplementsRFC7748(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.curve = X25519CurveGroupLaw()

    def test_rfc7748_iterative(self):
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
            1_000: "684cf59ba83309552800ef566f2f4d3c1c3887c49360e3875f2eb94d99532c51",
            1_000_000: "7c3911e0ab2586fd864497297e575e6f3bc601c0883c30df5f4dd2d24f665424",
        }

        for i in range(1, 1_000 + 1):
            k, u = crypto_scalarmult(k, u), k
            k_2, u_2 = self.curve.x25519(k_2, u_2), k_2
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
    
    def test_rfc7748_vectors(self):
        vectors = [(    
                "vector1",
                "a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4",
                "e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c",
                "c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552",
            ),
        ]

        for name, k_hex, u_hex, expected_hex in vectors:
            """Test X25519 against RFC 7748 test vectors."""
            k_bytes = unhexlify(k_hex)
            u_bytes = unhexlify(u_hex)
            expected = unhexlify(expected_hex)

            try:
                result_custom = self.curve.x25519(k_bytes, u_bytes)
            except Exception as e:
                self.fail(f"X25519 [{name}] implementation raised an exception: {e}")

            self.assertEqual(
                result_custom,
                expected,
                f"{name}: X25519 [{name}] implementation does not match the expected RFC 7748 result.",
            )

            self.assertEqual(
                result_custom, crypto_scalarmult(k_bytes, u_bytes),
                f"{name}: X25519 [{name}] implementation does not match PyNaCl's implementation."
            )