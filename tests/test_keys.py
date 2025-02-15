import secrets
import unittest

from parameterized import parameterized

from keys import KEY_SIZE, Key, PrivateKey, PublicKey, SharedKey


class TestKeys(unittest.TestCase):
    key_classes = [
        ("Key", Key),
        ("PrivateKey", PrivateKey),
        ("PublicKey", PublicKey),
        ("SharedKey", SharedKey),
    ]

    @parameterized.expand(key_classes)  # type: ignore
    def test_generated_key(self, name: str, key_class: type[Key]) -> None:
        """Test that a key is automatically generated and has the correct length."""
        instance = key_class()  # No key provided, so a random key is generated.
        self.assertIsInstance(instance.get_key(), bytes)
        self.assertEqual(
            len(instance.get_key()), KEY_SIZE, f"Key {name} is not the correct length."
        )

    @parameterized.expand(key_classes)  # type: ignore
    def test_manual_valid_key(self, name: str, key_class: type[Key]) -> None:
        """Test that providing a valid key works as expected."""
        valid_key = secrets.token_bytes(KEY_SIZE)
        instance = key_class(valid_key)
        self.assertEqual(
            instance.get_key(), valid_key, f"Key {name} is not equal to the provided key."
        )

    @parameterized.expand(key_classes)  # type: ignore
    def test_manual_invalid_key(self, name: str, key_class: type[Key]) -> None:
        """Test that providing an invalid key raises a ValueError."""
        invalid_key = secrets.token_bytes(KEY_SIZE - 1)
        with self.assertRaises(
            ValueError,
            msg=f"Key {name} should raise a ValueError when provided with invalid key.",
        ):
            _ = key_class(invalid_key)


if __name__ == "__main__":
    unittest.main()
