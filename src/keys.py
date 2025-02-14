import secrets
from dataclasses import dataclass

KEY_SIZE = 32


@dataclass
class Key:
    _key: bytes

    def __init__(self, key: bytes | None = None):
        self._key = secrets.token_bytes(KEY_SIZE) if key is None else key

    def __post_init__(self) -> None:
        if len(self._key) != KEY_SIZE:
            raise ValueError("Key must be 32 bytes long")

    def get_key(self) -> bytes:
        return self._key


@dataclass
class PrivateKey(Key):
    def __init__(self, key: bytes | None = None):
        super().__init__(key)


@dataclass
class PublicKey(Key): ...


@dataclass
class SharedKey(Key): ...
