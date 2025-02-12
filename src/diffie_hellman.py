from abc import ABC, abstractmethod
import os
from typing import Optional


class DiffieHellman(ABC):
    def __init__(self, private_key: Optional[bytes] = None):
        self.private_key = private_key if private_key is not None else os.urandom(32)

    @abstractmethod
    def generate_shared_secret(self, peer_public_key: bytes) -> bytes:
        raise NotImplementedError
    
    @abstractmethod
    def compute_public_key(self) -> bytes:
        raise NotImplementedError