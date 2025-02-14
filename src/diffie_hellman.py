from abc import ABC, abstractmethod

from keys import PrivateKey, PublicKey, SharedKey


class DiffieHellman(ABC):
    def __init__(self, private_key: PrivateKey):
        self.private_key = private_key

    @abstractmethod
    def generate_shared_secret(self, peer_public_key: PublicKey) -> SharedKey:
        raise NotImplementedError

    @abstractmethod
    def _compute_public_key(self) -> PublicKey:
        raise NotImplementedError
