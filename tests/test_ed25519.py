import secrets
import unittest
from parameterized import parameterized
from ed25519.affine_edwards_curve import AffineEdwardsCurve
from ed25519.edwards_curve import EdwardsCurve
from ed25519.extended_edwards_curve import ExtendedEdwardsCurve
from nacl.exceptions import BadSignatureError
from nacl.signing import SigningKey

from ed25519.edwards_signature_scheme import Ed25519


class TestEd25519Implementation(unittest.TestCase):
    
    def generate_keys(self, curve: EdwardsCurve) -> tuple[bytes, SigningKey, Ed25519]:
        """
        Generate a random 32-byte seed and create both a PyNaCl SigningKey and an instance
        of our custom Ed25519 implementation.

        Also, check that both implementations yield the same public key.
        """

        seed = secrets.token_bytes(32)
        nacl_signing_key = SigningKey(seed)
        nacl_public_key = nacl_signing_key.verify_key.encode()  # 32-byte public key

        custom_signer = Ed25519(secret_key=seed)

        # Check that both implementations produce the same public key.
        self.assertEqual(
            custom_signer.get_public_key(),
            nacl_public_key,
            "Public keys do not match between implementations.",
        )
        return seed, nacl_signing_key, custom_signer

    @parameterized.expand([
        (AffineEdwardsCurve(), "AffineEdwardsCurve"),
        (ExtendedEdwardsCurve(), "ExtendedEdwardsCurve"),
    ])
    def test_known_message(self, curve: EdwardsCurve, curve_name: str) -> None:
        _, nacl_signing_key, custom_signer = self.generate_keys(curve)
        msg = b"Attack at Dawn"

        nacl_signed = nacl_signing_key.sign(msg)
        nacl_signature = nacl_signed.signature  # 64-byte signature

        # Sign with our implementation
        custom_signature = custom_signer.sign(msg)

        # Signatures should match
        self.assertEqual(
            nacl_signature,
            custom_signature,
            f"Signatures do not match for a known message with {curve_name}.",
        )

        # Verify with PyNaCl (expects concat of signature + message)
        recovered_msg = nacl_signing_key.verify_key.verify(nacl_signature + msg)
        self.assertEqual(recovered_msg, msg)

        self.assertTrue(
            custom_signer.verify(custom_signature, msg, custom_signer.public_key),
            f"Custom verification failed for a known message with {curve_name}.",
        )
        
    @parameterized.expand([
        (AffineEdwardsCurve(), "AffineEdwardsCurve"),
        (ExtendedEdwardsCurve(), "ExtendedEdwardsCurve"),
    ])
    def test_empty_message(self, curve: EdwardsCurve, curve_name: str) -> None:
        # Test signing and verifying an empty message.
        _, nacl_signing_key, custom_signer = self.generate_keys(curve)
        msg = b""

        nacl_signature = nacl_signing_key.sign(msg).signature
        custom_signature = custom_signer.sign(msg)

        self.assertEqual(
            nacl_signature,
            custom_signature,
            f"Signatures do not match for an empty message with {curve_name}.",
        )

        # Verify using both implementations
        recovered_msg = nacl_signing_key.verify_key.verify(nacl_signature + msg)
        self.assertEqual(recovered_msg, msg)
        self.assertTrue(
            custom_signer.verify(custom_signature, msg, custom_signer.public_key),
            f"Custom verification failed for an empty message with {curve_name}.",
        )

    @parameterized.expand([
        (AffineEdwardsCurve(), "AffineEdwardsCurve"),
        (ExtendedEdwardsCurve(), "ExtendedEdwardsCurve"),
    ])
    def test_random_messages(self, curve: EdwardsCurve, curve_name: str) -> None:
        # Generate and test a variety of random messages.
        _, nacl_signing_key, custom_signer = self.generate_keys(curve)
        for _ in range(100):
            # Create a message of random length between 0 and 1024 bytes.
            msg = secrets.token_bytes(secrets.randbelow(1025))
            nacl_signature = nacl_signing_key.sign(msg).signature
            custom_signature = custom_signer.sign(msg)

            self.assertEqual(
                nacl_signature,
                custom_signature,
                f"Signatures do not match for message: {msg.hex()} with {curve_name}.",
            )

            recovered_msg = nacl_signing_key.verify_key.verify(nacl_signature + msg)
            self.assertEqual(recovered_msg, msg)
            self.assertTrue(
                custom_signer.verify(custom_signature, msg, custom_signer.public_key),
                f"Custom verification failed for message: {msg.hex()}",
            )

    @parameterized.expand([
        (AffineEdwardsCurve(), "AffineEdwardsCurve"),
        (ExtendedEdwardsCurve(), "ExtendedEdwardsCurve"),
    ])
    def test_invalid_signature(self, curve: EdwardsCurve, curve_name: str) -> None:
        _, nacl_signing_key, custom_signer = self.generate_keys(curve)
        msg = b"Test message for invalid signature"
        custom_signature = custom_signer.sign(msg)

        altered_signature = bytearray(custom_signature)
        altered_signature[0] ^= 0x01  # flip a bit in the first byte

        # PyNaCl verification should raise a BadSignatureError
        with self.assertRaises(BadSignatureError):
            nacl_signing_key.verify_key.verify(bytes(altered_signature) + msg)

        self.assertFalse(
            custom_signer.verify(altered_signature, msg, custom_signer.public_key),
            f"Custom verification incorrectly accepted an altered signature with {curve_name}.",
        )

    @parameterized.expand([
        (AffineEdwardsCurve(), "AffineEdwardsCurve"),
        (ExtendedEdwardsCurve(), "ExtendedEdwardsCurve"),
    ])
    def test_wrong_public_key(self, curve: EdwardsCurve, curve_name: str) -> None:
        _, _, custom_signer = self.generate_keys(curve)
        msg = b"Test message with wrong public key"
        custom_signature = custom_signer.sign(msg)

        wrong_seed = secrets.token_bytes(32)
        wrong_signer = Ed25519(secret_key=wrong_seed)

        self.assertNotEqual(
            custom_signer.public_key,
            wrong_signer.public_key,
            f"Unexpectedly, both signers have the same public key with {curve_name}.",
        )

        self.assertFalse(
            custom_signer.verify(custom_signature, msg, wrong_signer.public_key),
            f"Custom verify accepted a signature with the wrong public key with {curve_name}.",
        )

    @parameterized.expand([
        (AffineEdwardsCurve(), "AffineEdwardsCurve"),
        (ExtendedEdwardsCurve(), "ExtendedEdwardsCurve"),
    ])
    def test_large_message(self, curve: EdwardsCurve, curve_name: str) -> None:
        # Test with a large message (e.g. 10 KB).
        _, nacl_signing_key, custom_signer = self.generate_keys(curve)
        msg = secrets.token_bytes(10 * 1024)

        nacl_signature = nacl_signing_key.sign(msg).signature
        custom_signature = custom_signer.sign(msg)

        self.assertEqual(
            nacl_signature,
            custom_signature,
            f"Signatures do not match for a large message with {curve_name}.",
        )

        recovered_msg = nacl_signing_key.verify_key.verify(nacl_signature + msg)
        self.assertEqual(recovered_msg, msg)
        self.assertTrue(
            custom_signer.verify(custom_signature, msg, custom_signer.public_key),
            f"Custom verification failed for a large message with {curve_name}.",
        )

    
    # TODO: CS run timing tests

if __name__ == "__main__":
    unittest.main()