import hashlib

from ed25519.edwards_curve import EdwardsCurve
from util import clamp_scalar


class Ed25519:
    """
    An implementation of the Ed25519 signature scheme.

    This class uses an EdwardsCurve instance for all curve arithmetic.
    """

    def __init__(self, secret_key: bytes):
        self.curve = EdwardsCurve()
        self._hashed_secret_key = hashlib.sha512(secret_key).digest()
        s_bits = bytearray(self._hashed_secret_key[:32])
        s_int = clamp_scalar(s_bits)
        pk_point = self.curve.scalar_mult(s_int, self.curve.B)
        self.public_key = self.curve.compress(pk_point)

    def sign(self, msg: bytes) -> bytes:
        """
        Sign a message using Ed25519.

        Steps:
          1. Compute h = SHA512(secret_key) → (s_bits, prefix) (each 32 bytes).
          2. Clamp s_bits to get the scalar s.
          3. Compute public key: pk = compress(s * B).
          4. Compute nonce r = SHA512(prefix || msg) mod q.
          5. Compute R = compress(r * B).
          6. Compute challenge k = SHA512(R || pk || msg) mod q.
          7. Compute response t = (r + k * s) mod q.
          8. Return signature = R || t (64 bytes).
        """
        # Split the hashed secret key into its two halves.
        s_bits = bytearray(self._hashed_secret_key[:32])
        prefix = self._hashed_secret_key[32:]
        s_int = clamp_scalar(s_bits)

        # Compute nonce r = SHA512(prefix || msg) mod q.
        r_hash = hashlib.sha512(prefix + msg).digest()
        r_int = int.from_bytes(r_hash, "little") % self.curve.q
        R_point = self.curve.scalar_mult(r_int, self.curve.B)
        R_comp = self.curve.compress(R_point)

        # Compute challenge k = SHA512(R || public_key || msg) mod q.
        k_hash = hashlib.sha512(R_comp + self.public_key + msg).digest()
        k_int = int.from_bytes(k_hash, "little") % self.curve.q

        # Compute response t = (r + k * s) mod q.
        t_int = (r_int + k_int * s_int) % self.curve.q
        t_bytes = t_int.to_bytes(32, "little")

        return R_comp + t_bytes

    def verify(self, sig: bytes, msg: bytes, pk: bytes) -> bool:
        """
        Verify an Ed25519 signature.

        Steps:
          1. Split sig into R (first 32 bytes) and t (last 32 bytes).
          2. Convert t to an integer modulo q.
          3. Uncompress R and pk to get the corresponding curve points.
          4. Compute challenge k = SHA512(R || pk || msg) mod q.
          5. Check that [t]B equals R + [k]A (where A is the public key point).
        """
        curve = self.curve

        if len(sig) != 64:
            raise ValueError("Signature must be 64 bytes")
        if len(pk) != 32:
            raise ValueError("Public key must be 32 bytes")
        R_comp = sig[:32]
        t_bytes = sig[32:]
        t_int = int.from_bytes(t_bytes, "little") % curve.q

        try:
            R = curve.uncompress(R_comp)
            A = curve.uncompress(pk)
        except ValueError:
            return False

        k_hash = hashlib.sha512(R_comp + pk + msg).digest()
        k_int = int.from_bytes(k_hash, "little") % curve.q

        # Compute left-hand side: [t]B.
        LHS = curve.scalar_mult(t_int, curve.B)
        # Compute right-hand side: R + [k]A.
        kA = curve.scalar_mult(k_int, A)
        RHS = curve.add(R, kA)

        return LHS == RHS
