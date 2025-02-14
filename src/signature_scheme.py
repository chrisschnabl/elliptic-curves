from abc import ABC, abstractmethod

from keys import PrivateKey, PublicKey


class SignatureScheme(ABC):
    secret_key: PrivateKey
    public_key: PublicKey

    @abstractmethod
    def sign(self, msg: bytes) -> bytes:
        pass

    @abstractmethod
    def verify(self, sig: bytes, msg: bytes, pk: PublicKey) -> bool:
        pass

    @abstractmethod
    def get_public_key(self) -> bytes:
        pass
