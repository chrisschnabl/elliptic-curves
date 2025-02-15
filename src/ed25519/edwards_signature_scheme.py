from collections.abc import Callable

import nacl.hash

from ed25519.affine_edwards_curve import AffineEdwardsCurve
from ed25519.edwards_curve import EdwardsCurve
from keys import PrivateKey, PublicKey
from signature_scheme import SignatureScheme
from util import clamp_scalar


class Ed25519(SignatureScheme):  # type: ignore
    """
    An implementation of the Ed25519 signature scheme.

    This class uses an EdwardsCurve instance for all curve arithmetic.
    """

    def __init__(
        self, secret_key: PrivateKey, curve: EdwardsCurve = AffineEdwardsCurve()
    ):
        self.curve: EdwardsCurve = curve

        def hash_function(plain_text: bytes) -> bytes:
            return nacl.hash.sha512(plain_text, encoder=nacl.encoding.RawEncoder)  # type: ignore

        self.hash_function: Callable[[bytes], bytes] = hash_function
        self._hashed_secret_key = self.hash_function(secret_key.get_key())
        s_bits = self._hashed_secret_key[:32]
        self.s_int = clamp_scalar(bytearray(s_bits))
        self.public_key = self.curve.scalar_mult(self.curve.B, self.s_int)
        self.public_key = PublicKey(self.curve.compress(self.public_key))

    def sign(self, msg: bytes) -> bytes:
        """
        Sign a message using Ed25519.

        Steps (taken from the Slide):
          1. Compute h = SHA512(secret_key) -> (s_bits, prefix) (each 32 bytes).
          2. Clamp s_bits to get the scalar s.
          3. Compute and compress public key: pk = compress(s * B).
          4. Compute r = SHA512(prefix || msg) mod q. (q is the order of the Base Point)
          5. Compute R = compress(r * B).
          6. Compute challenge k = SHA512(R || pk || msg) mod q.
          7. Compute response t = (r + k * s) mod q.
          8. Return signature = R || t (64 bytes).
        """
        prefix = self._hashed_secret_key[32:]

        # Compute r = SHA512(prefix || msg) mod q.
        r_hash = self.hash_function(prefix + msg)
        r_int = int.from_bytes(r_hash, "little") % self.curve.q
        R_point = self.curve.scalar_mult(self.curve.B, r_int)
        R_comp = self.curve.compress(R_point)

        # Compute challenge k = SHA512(R || public_key || msg) mod q.
        k_hash = self.hash_function(R_comp + self.public_key.get_key() + msg)
        k_int = int.from_bytes(k_hash, "little") % self.curve.q

        # Compute response t = (r + k * s) mod q.
        t_int = (r_int + k_int * self.s_int) % self.curve.q
        t_bytes = t_int.to_bytes(32, "little")

        return R_comp + t_bytes  # type: ignore

    def verify(self, sig: bytes, msg: bytes, pk: PublicKey) -> bool:
        """
        Verify an Ed25519 signature.

        Steps:
          1. Split sig into R (first 32 bytes) and t (last 32 bytes).
          2. Convert t to an integer modulo q.
          3. Uncompress R and pk to get the corresponding curve points.
          4. Compute challenge k = SHA512(R || pk || msg) mod q.
          5. Check that [t]B equals R + [k]A (where A is the public key point).
        """

        if len(sig) != 64:
            raise ValueError("Signature must be 64 bytes")

        R_comp = sig[:32]
        t_bytes = sig[32:]
        t_int = int.from_bytes(t_bytes, "little") % self.curve.q

        try:
            R = self.curve.uncompress(R_comp)
            A = self.curve.uncompress(pk.get_key())
        except ValueError:
            return False

        k_hash = self.hash_function(R_comp + pk.get_key() + msg)
        k_int = int.from_bytes(k_hash, "little") % self.curve.q

        # Compute left-hand side: [t]B.
        LHS = self.curve.scalar_mult(
            self.curve.B, t_int
        )  # B is extended homogeneous coordinates
        # Compute right-hand side: R + [k]A.
        kA = self.curve.scalar_mult(A, k_int)
        RHS = self.curve.add(R, kA)

        return self.curve.point_equals(LHS, RHS)  # type: ignore

    def get_public_key(self) -> PublicKey:
        return self.public_key
