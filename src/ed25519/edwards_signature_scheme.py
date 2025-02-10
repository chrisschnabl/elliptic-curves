from ed25519.affine_edwards_curve import AffineEdwardsCurve
from ed25519.edwards_curve import EdwardsCurve
from ed25519.extended_edwards_curve import ExtendedEdwardsCurve
from ed25519.signature_scheme import SignatureScheme
import nacl.hash
from util import clamp_scalar
import hashlib
class Ed25519(SignatureScheme):
    """
    An implementation of the Ed25519 signature scheme.

    This class uses an EdwardsCurve instance for all curve arithmetic.
    """

    def __init__(self, secret_key: bytes, curve: EdwardsCurve = AffineEdwardsCurve()):
        self.curve = curve
        self.hash_function = lambda plain_text: hashlib.sha512(plain_text).digest()
        #self.hash_function = lambda plain_text: nacl.hash.sha512(plain_text, encoder=nacl.encoding.RawEncoder)
        # somehow this fails for the invalid Signature TEST TODO CS: why?
        self._hashed_secret_key = self.hash_function(secret_key)
        self.s_int = clamp_scalar(bytearray(self._hashed_secret_key[:32]))
        self.public_key = self.curve.scalar_mult(self.curve.B, self.s_int)
        self.public_key = self.curve.compress(self.public_key)


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
        prefix = self._hashed_secret_key[32:]

        # Compute nonce r = SHA512(prefix || msg) mod q.
        r_hash = self.hash_function(prefix + msg)
        r_int = int.from_bytes(r_hash, "little") % self.curve.q
        R_point = self.curve.scalar_mult(self.curve.B, r_int)
        R_comp = self.curve.compress(R_point)

        # Compute challenge k = SHA512(R || public_key || msg) mod q.
        k_hash = self.hash_function(R_comp + self.public_key + msg)
        k_int = int.from_bytes(k_hash, "little") % self.curve.q

        # Compute response t = (r + k * s) mod q.
        t_int = (r_int + k_int * self.s_int) % self.curve.q
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

        if len(sig) != 64:
            raise ValueError("Signature must be 64 bytes")
        if len(pk) != 32:
            raise ValueError("Public key must be 32 bytes")
        R_comp = sig[:32]
        t_bytes = sig[32:]
        t_int = int.from_bytes(t_bytes, "little") % self.curve.q

        try:
            R = self.curve.uncompress(R_comp)
            A = self.curve.uncompress(pk)
        except ValueError:
            return False

        k_hash = self.hash_function(R_comp + pk + msg)
        k_int = int.from_bytes(k_hash, "little") % self.curve.q

        # Compute left-hand side: [t]B.
        LHS = self.curve.scalar_mult(self.curve.B, t_int) # B is extended homogeneous coordinates
        # Compute right-hand side: R + [k]A.
        kA = self.curve.scalar_mult(A, k_int)
        RHS = self.curve.add(R, kA)

        return self.curve.point_equals(LHS, RHS)  # Fails if lambda scalar scales the projective coordinates
    
    def get_public_key(self) -> bytes:
        return self.public_key