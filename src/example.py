from ed25519.edwards_signature_scheme import Ed25519
from keys import PrivateKey, PublicKey
from x25519.elliptic_curve_diffie_hellman import EllipticCurveDiffieHellman

# This is just a simple example of how to use the library.
# This is not meant to be production ready code.


def minimal_dh_example() -> None:
    my_private_key = PrivateKey()
    dh = EllipticCurveDiffieHellman(my_private_key)
    my_public_key = dh.compute_public_key()

    peer_private_key = PrivateKey()
    dh_peer = EllipticCurveDiffieHellman(peer_private_key)
    peer_public_key = dh_peer.compute_public_key()
    dh.generate_shared_secret(peer_public_key)
    dh_peer.generate_shared_secret(my_public_key)


def minimal_ed25519_example() -> None:
    my_private_key = PrivateKey()
    signer = Ed25519(secret_key=my_private_key)

    message = b"Hello, Ed25519!"
    signer.sign(message)
    pub_key_bytes = signer.get_public_key()
    PublicKey(pub_key_bytes.get_key())


if __name__ == "__main__":
    minimal_dh_example()
    minimal_ed25519_example()
