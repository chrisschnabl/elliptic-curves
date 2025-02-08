import secrets
import unittest

from nacl.exceptions import BadSignatureError
from nacl.signing import SigningKey

from ed25519.ed25519_signatures import Ed25519


class TestEd25519Implementation(unittest.TestCase):
    def generate_keys(self) -> tuple[bytes, SigningKey, Ed25519]:
        """
        Generate a random 32-byte seed and create both a PyNaCl SigningKey and an instance
        of your custom Ed25519 implementation.

        Also, check that both implementations yield the same public key.
        """

        seed = secrets.token_bytes(32)
        nacl_signing_key = SigningKey(seed)
        nacl_public_key = nacl_signing_key.verify_key.encode()  # 32-byte public key

        # Create an instance of your custom signer.
        custom_signer = Ed25519(secret_key=seed)

        # Check that both implementations produce the same public key.
        self.assertEqual(
            custom_signer.public_key,
            nacl_public_key,
            "Public keys do not match between implementations.",
        )
        return seed, nacl_signing_key, custom_signer

    def test_known_message(self) -> None:
        # Test with a fixed message ('Attack at Dawn').
        _, nacl_signing_key, custom_signer = self.generate_keys()
        msg = b"Attack at Dawn"

        # Sign with PyNaCl
        nacl_signed = nacl_signing_key.sign(msg)
        nacl_signature = nacl_signed.signature  # 64-byte signature

        # Sign with our implementation
        custom_signature = custom_signer.sign(msg)

        # Signatures should match
        self.assertEqual(
            nacl_signature,
            custom_signature,
            "Signatures do not match for a known message.",
        )

        # Verify with PyNaCl (expects concat of signature + message)
        recovered_msg = nacl_signing_key.verify_key.verify(nacl_signature + msg)
        self.assertEqual(recovered_msg, msg)

        self.assertTrue(
            custom_signer.verify(custom_signature, msg, custom_signer.public_key),
            "Custom verification failed for a known message.",
        )

    def test_empty_message(self) -> None:
        # Test signing and verifying an empty message.
        _, nacl_signing_key, custom_signer = self.generate_keys()
        msg = b""

        nacl_signature = nacl_signing_key.sign(msg).signature
        custom_signature = custom_signer.sign(msg)

        self.assertEqual(
            nacl_signature,
            custom_signature,
            "Signatures do not match for an empty message.",
        )

        # Verify using both implementations
        recovered_msg = nacl_signing_key.verify_key.verify(nacl_signature + msg)
        self.assertEqual(recovered_msg, msg)
        self.assertTrue(
            custom_signer.verify(custom_signature, msg, custom_signer.public_key),
            "Custom verification failed for an empty message.",
        )

    def test_random_messages(self) -> None:
        # Generate and test a variety of random messages.
        _, nacl_signing_key, custom_signer = self.generate_keys()
        for _ in range(1):
            # Create a message of random length between 0 and 1024 bytes.
            msg = secrets.token_bytes(secrets.randbelow(1025))
            nacl_signature = nacl_signing_key.sign(msg).signature
            custom_signature = custom_signer.sign(msg)

            self.assertEqual(
                nacl_signature,
                custom_signature,
                f"Signatures do not match for message: {msg.hex()}",
            )

            recovered_msg = nacl_signing_key.verify_key.verify(nacl_signature + msg)
            self.assertEqual(recovered_msg, msg)
            self.assertTrue(
                custom_signer.verify(custom_signature, msg, custom_signer.public_key),
                f"Custom verification failed for message: {msg.hex()}",
            )

    def test_invalid_signature(self) -> None:
        # Verify that a modified (invalid) signature is rejected by both implementations.
        _, nacl_signing_key, custom_signer = self.generate_keys()
        msg = b"Test message for invalid signature"
        custom_signature = custom_signer.sign(msg)

        # Modify one byte in the signature (flip one bit).
        altered_signature = bytearray(custom_signature)
        altered_signature[0] ^= 0x01  # flip a bit in the first byte

        # PyNaCl verification should raise a BadSignatureError
        with self.assertRaises(BadSignatureError):
            nacl_signing_key.verify_key.verify(bytes(altered_signature) + msg)

        # Your custom implementation should return False.
        self.assertFalse(
            custom_signer.verify(altered_signature, msg, custom_signer.public_key),
            "Custom verification incorrectly accepted an altered signature.",
        )

    def test_wrong_public_key(self) -> None:
        # Test that using a wrong public key causes signature verification to fail.
        _, nacl_signing_key, custom_signer = self.generate_keys()
        msg = b"Test message with wrong public key"
        custom_signature = custom_signer.sign(msg)

        # Generate a different key pair.
        wrong_seed = secrets.token_bytes(32)
        wrong_signer = Ed25519(secret_key=wrong_seed)

        self.assertNotEqual(
            custom_signer.public_key,
            wrong_signer.public_key,
            "Unexpectedly, both signers have the same public key.",
        )

        # Verification using the wrong public key should fail.
        self.assertFalse(
            custom_signer.verify(custom_signature, msg, wrong_signer.public_key),
            "Custom verify accepted a signature with the wrong public key.",
        )

    def test_large_message(self) -> None:
        # Test with a large message (e.g. 10 KB).
        _, nacl_signing_key, custom_signer = self.generate_keys()
        msg = secrets.token_bytes(10 * 1024)

        nacl_signature = nacl_signing_key.sign(msg).signature
        custom_signature = custom_signer.sign(msg)

        self.assertEqual(
            nacl_signature,
            custom_signature,
            "Signatures do not match for a large message.",
        )

        recovered_msg = nacl_signing_key.verify_key.verify(nacl_signature + msg)
        self.assertEqual(recovered_msg, msg)
        self.assertTrue(
            custom_signer.verify(custom_signature, msg, custom_signer.public_key),
            "Custom verification failed for a large message.",
        )


if __name__ == "__main__":
    unittest.main()
