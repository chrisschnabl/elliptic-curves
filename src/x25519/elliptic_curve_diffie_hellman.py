from diffie_hellman import DiffieHellman
from keys import PrivateKey, PublicKey, SharedKey
from util import encode_u_coordinate
from x25519.curve25519 import Curve25519
from x25519.montgomery_ladder import MontgomeryLadderRFC7748


class EllipticCurveDiffieHellman(DiffieHellman):  # type: ignore
    """
    A Diffie=Hellman abstraction for the X25519 key exchange.

    This class uses a provided curve for scalar multiplication and:
      - Generates a random 32-byte private key (use better randomness here for production)
      - Derive the corresponding 32-byte public key.
      - Derive of a shared secret
    """

    def __init__(
        self, private_key: PrivateKey, curve: Curve25519 | None = None
    ):  # sub: ignore
        """
        Initialize the DiffieHellman instance.

        Args:
            private_key (Optional[bytes]): A 32-byte private key.
            ladder (Optional[MontgomeryLadder]): An instance of the Montgomery ladder.
        """
        super().__init__(private_key)
        self.curve = curve if curve is not None else MontgomeryLadderRFC7748()
        self.public_key = self.compute_public_key()

    def compute_public_key(self) -> PublicKey:
        """
        Computes the public key corresponding to the private key.

        For X25519 the standard base point is defined as 9, represented as 0x09
        followed by 31 zero bytes (little-endian).

        Returns:
            PublicKey: The 32-byte public key.
        """
        return PublicKey(
            self.curve.x25519(
                self.private_key.get_key(), encode_u_coordinate(self.curve.B.x)
            )
        )

    def generate_shared_secret(self, peer_public_key: PublicKey) -> SharedKey:
        """
        Computes the shared secret given a peer's public key.

        Args:
            peer_public_key (PublicKey): The peer's 32-byte public key.

        Returns:
            SharedKey: The computed 32-byte shared secret.
        """
        return SharedKey(
            self.curve.x25519(self.private_key.get_key(), peer_public_key.get_key())
        )
