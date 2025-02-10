from abc import ABC, abstractmethod

class SignatureScheme(ABC):
    
    def __init__(self):
        self.secret_key = None
        
    @abstractmethod
    def sign(self, msg: bytes) -> bytes:
        pass

    @abstractmethod
    def verify(self, sig: bytes, msg: bytes, pk: bytes) -> bool:
        pass

    @abstractmethod
    def get_public_key(self) -> bytes:
        pass