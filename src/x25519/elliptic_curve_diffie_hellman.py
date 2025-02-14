from diffie_hellman import DiffieHellman
from keys import PrivateKey, PublicKey, SharedKey
from x25519.montgomery_ladder import MontgomeryLadderRFC7748
from x25519.x25519_curve import X25519Curve


class EllipticCurveDiffieHellman(DiffieHellman):  # type: ignore
    """
    A Diffie=Hellman abstraction for the X25519 key exchange.

    This class uses a provided curve for scalar multiplication and:
      - Generates a random 32-byte private key (use better randomness here for production)
      - Derive the corresponding 32-byte public key.
      - Derive of a shared secret
    """

    def __init__(
        self, private_key: PrivateKey, curve: X25519Curve | None = None
    ):  # sub: ignore
        """
        Initialize the DiffieHellman instance.

        Args:
            private_key (Optional[bytes]): A 32-byte private key.
            ladder (Optional[MontgomeryLadder]): An instance of the Montgomery ladder.
        """
        super().__init__(private_key)
        self.x25519 = curve if curve is not None else MontgomeryLadderRFC7748()
        self.public_key = self._compute_public_key()

    def _compute_public_key(self) -> PublicKey:
        """
        Computes the public key corresponding to the private key.

        For X25519 the standard base point is defined as 9, represented as 0x09
        followed by 31 zero bytes (little-endian).

        Returns:
            PublicKey: The 32-byte public key.
        """
        base_point = b"\x09" + (b"\x00" * 31)  # TODO CS: re-use the curve's base point
        return PublicKey(self.x25519.x25519(self.private_key.get_key(), base_point))

    def generate_shared_secret(self, peer_public_key: PublicKey) -> SharedKey:
        """
        Computes the shared secret given a peer's public key.

        Args:
            peer_public_key (PublicKey): The peer's 32-byte public key.

        Returns:
            SharedKey: The computed 32-byte shared secret.
        """
        return SharedKey(
            self.x25519.x25519(self.private_key.get_key(), peer_public_key.get_key())
        )


# -------------------------------------------------------------------
# Example Usage
# -------------------------------------------------------------------

if __name__ == "__main__":
    # Party A: Alice creates her Diffie–Hellman object (with a random private key)
    alice = EllipticCurveDiffieHellman(PrivateKey())

    # Party B: Bob creates his Diffie–Hellman object
    bob = EllipticCurveDiffieHellman(PrivateKey())

    # Each party computes the shared secret using the other's public key.
    alice_shared = alice.generate_shared_secret(bob.public_key)
    bob_shared = bob.generate_shared_secret(alice.public_key)

    # The shared secrets should match.
    match = alice_shared == bob_shared
    if not match:
        pass
    else:
        pass
