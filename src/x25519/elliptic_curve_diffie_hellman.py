from typing import Optional
from diffie_hellman import DiffieHellman
from x25519.montgomery_ladder import MontgomeryLadderRFC7748
from x25519.x25519_curve import X25519Curve


class EllipticCurveDiffieHellman(DiffieHellman):
    """
    A Diffie=Hellman abstraction for the X25519 key exchange.

    This class uses a Montgomery ladder for scalar multiplication to provide:
      - Generate a random 32-byte private key (use better randomness here for production)
      - Derive the corresponding 32-byte public key.
      - Derive of a shared secret
    """

    def __init__(self, private_key: Optional[bytes] = None,
                 curve: Optional[X25519Curve] = None):
        """
        Initialize the DiffieHellman instance.

        Args:
            private_key (Optional[bytes]): A 32-byte private key. If not provided,
                                           one will be generated using os.urandom.
            ladder (Optional[MontgomeryLadder]): An instance of the Montgomery ladder.
        """
        super().__init__(private_key)
        self.x25519 = curve if curve is not None else MontgomeryLadderRFC7748()
        self.public_key = self.compute_public_key()

    def compute_public_key(self) -> bytes:
        """
        Computes the public key corresponding to the private key.

        For X25519 the standard base point is defined as 9, represented as 0x09
        followed by 31 zero bytes (little-endian).

        Returns:
            bytes: The 32-byte public key.
        """
        base_point = b'\x09' + (b'\x00' * 31)
        return self.x25519.x25519(self.private_key, base_point)

    def generate_shared_secret(self, peer_public_key: bytes) -> bytes:
        """
        Computes the shared secret given a peer's public key.

        Args:
            peer_public_key (bytes): The peer's 32-byte public key.

        Returns:
            bytes: The computed 32-byte shared secret.
        """
        return self.x25519.x25519(self.private_key, peer_public_key)
    

# -------------------------------------------------------------------
# Example Usage
# -------------------------------------------------------------------

if __name__ == "__main__":
    # Party A: Alice creates her Diffie–Hellman object (with a random private key)
    alice = EllipticCurveDiffieHellman()
    print("Alice's Public Key:", alice.public_key.hex())

    # Party B: Bob creates his Diffie–Hellman object
    bob = EllipticCurveDiffieHellman()
    print("Bob's Public Key:", bob.public_key.hex())

    # Each party computes the shared secret using the other's public key.
    alice_shared = alice.generate_shared_secret(bob.public_key)
    bob_shared = bob.generate_shared_secret(alice.public_key)

    print("Alice's Shared Secret:", alice_shared.hex())
    print("Bob's Shared Secret:", bob_shared.hex())

    # The shared secrets should match.
    assert alice_shared == bob_shared, "Shared secrets do not match!"