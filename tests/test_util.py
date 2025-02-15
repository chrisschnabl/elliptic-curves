import secrets
import unittest

from util import (
    affine_to_projective,
    clamp_scalar,
    cswap,
    decode_u,
    encode_u_coordinate,
    modinv,
    projective_to_affine,
    sqrt_mod,
)

# Import the functions from your module.
# For example, if they are defined in a module called `crypto_math`:
# from crypto_math import (
#     cswap, clamp_scalar, modinv, sqrt_mod,
#     decode_u, encode_u_coordinate,
#     affine_to_projective, projective_to_affine
# )

# For this example, we'll assume the functions are in the current namespace.


# --- Tests for cswap ---
class TestCSwap(unittest.TestCase):
    def test_no_swap(self) -> None:
        p = 101  # A small prime
        x = 37
        y = 73
        swap = 0
        x_new, y_new = cswap(swap, x, y, p)
        self.assertEqual(x_new, x % p)
        self.assertEqual(y_new, y % p)

    def test_swap(self) -> None:
        p = 101
        x = 37
        y = 73
        swap = 1
        x_new, y_new = cswap(swap, x, y, p)
        self.assertEqual(x_new, y % p)
        self.assertEqual(y_new, x % p)

    def test_modulo_behavior(self) -> None:
        # Test with values larger than p to ensure the mod is applied.
        p = 97
        x = 150  # 150 mod 97 = 53
        y = 80  # 80 mod 97 = 80
        swap = 1
        x_new, y_new = cswap(swap, x, y, p)
        self.assertEqual(x_new, (y % p))
        self.assertEqual(y_new, (x % p))


# --- Tests for clamp_scalar ---
class TestClampScalar(unittest.TestCase):
    def test_clamp_scalar_all_ones(self) -> None:
        # Input: 32 bytes of 0xff.
        k = bytearray(b"\xff" * 32)
        result = clamp_scalar(k)
        # After clamping:
        # - k[0] becomes 0xff & 248 = 0xf8.
        self.assertEqual(k[0], 0xF8)
        # - k[31] becomes (0xff & 127) | 64 = (127 | 64) = 127 (0x7F).
        self.assertEqual(k[31], 0x7F)
        # And the returned integer should equal int.from_bytes(k, "little")
        self.assertEqual(result, int.from_bytes(k, "little"))

    def test_clamp_scalar_all_zero(self) -> None:
        k = bytearray(b"\x00" * 32)
        result = clamp_scalar(k)
        # k[0] remains 0; k[31] becomes (0 & 127) | 64 = 64.
        self.assertEqual(k[0], 0x00)
        self.assertEqual(k[31], 64)
        self.assertEqual(result, int.from_bytes(k, "little"))


# --- Tests for modinv ---
class TestModInv(unittest.TestCase):
    def test_modinv_nonzero(self) -> None:
        p = 101
        # For each nonzero x mod p, check that (x * modinv(x, p)) % p == 1.
        for x in range(1, p):
            inv = modinv(x, p)
            self.assertEqual((x * inv) % p, 1)

    def test_modinv_zero(self) -> None:
        # By our implementation, modinv(0, p) returns 0,
        # even though mathematically 0 has no inverse.
        p = 101
        self.assertEqual(modinv(0, p), 0)


# --- Tests for sqrt_mod ---
class TestSqrtMod(unittest.TestCase):
    def test_sqrt_mod_zero(self) -> None:
        p = 13
        self.assertEqual(sqrt_mod(0, p), 0)

    def test_sqrt_mod_branch1(self) -> None:
        # For some quadratic residue that satisfies r^2 ≡ a mod p directly.
        p = 13
        a = 3  # 3 is a quadratic residue modulo 13 (since 4^2=16≡3 mod 13, for example)
        r = sqrt_mod(a, p)
        self.assertEqual((r * r) % p, a % p)

    def test_sqrt_mod_branch2(self) -> None:
        # Choose a value that triggers the second branch.
        # For p = 13, let a = 4.
        # Then, exp = (13+3)//8 = 2, and r = 4**2 mod 13 = 16 mod 13 = 3.
        # Since 3^2 = 9 and 9 ≡ -4 mod 13, branch two should be taken.
        # sqrt_m1 = 2^((13-1)//4) mod 13 = 2^3 mod 13 = 8.
        # Expected r2 = (3 * 8) mod 13 = 24 mod 13 = 11.
        p = 13
        a = 4
        r = sqrt_mod(a, p)
        self.assertEqual((r * r) % p, a % p)
        # r should be either 11 or its negative modulo 13 (i.e. 2)
        self.assertIn(r, [11, 2])

    def test_sqrt_mod_no_root(self) -> None:
        # For p = 13, choose a value that is not a quadratic residue.
        p = 13
        a = 2  # 2 is not a quadratic residue mod 13.
        with self.assertRaises(ValueError):
            sqrt_mod(a, p)


# --- Tests for decode_u and encode_u_coordinate ---
class TestDecodeEncodeU(unittest.TestCase):
    def test_decode_u(self) -> None:
        # Create a 32-byte string with nonzero values and with
        # the top bit of the last byte set.
        u_bytes = bytearray(b"\x01" * 31 + b"\xff")
        decoded = decode_u(bytes(u_bytes))
        # The function masks off the high bit of the last byte.
        expected_bytes = bytearray(u_bytes)
        expected_bytes[31] &= 127
        expected_int = int.from_bytes(expected_bytes, "little")
        self.assertEqual(decoded, expected_int)

    def test_encode_u_coordinate(self) -> None:
        x = 12345678901234567890
        encoded = encode_u_coordinate(x)
        self.assertEqual(len(encoded), 32)
        # Decoding the encoded coordinate should yield x.
        self.assertEqual(int.from_bytes(encoded, "little"), x)

    def test_round_trip(self) -> None:
        # Test that encoding an integer and converting it back produces the same integer.
        x = secrets.randbits(256)
        encoded = encode_u_coordinate(x)
        decoded = int.from_bytes(encoded, "little")
        self.assertEqual(x, decoded)


# --- Tests for affine/projective conversions ---
class TestAffineProjective(unittest.TestCase):
    def test_affine_to_projective(self) -> None:
        x = 42
        X, Z = affine_to_projective(x)
        self.assertEqual(X, x)
        self.assertEqual(Z, 1)

    def test_projective_to_affine(self) -> None:
        p = 101
        # For a given x and nonzero Z, if we set X = x * Z mod p,
        # then converting back should yield x mod p.
        for x in range(0, p):
            for Z in range(1, p):  # Z must be invertible mod p.
                X = (x * Z) % p
                x_affine = projective_to_affine(X, Z, p)
                self.assertEqual(x_affine, x % p)


if __name__ == "__main__":
    unittest.main()
