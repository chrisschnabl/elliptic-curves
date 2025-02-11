from abc import ABC, abstractmethod
from typing import Callable

from util import clamp_scalar, decode_u, encode_u_coordinate

class X25519Implementation(ABC):
    def __init__(self):
       self.p = 2**255 - 19
       self.a24 = 121665

    def __call__(self, k_bytes: bytes, u_bytes: bytes) -> bytes:
        return self.ladder(k_bytes, u_bytes)  
       
    @abstractmethod
    def ladder(self, k_int: int, xP: int) -> int:
        raise NotImplementedError

class EllipticCurveDiffieHellman:
    def __init__(self, ladder: X25519Implementation):
      self.ladder = ladder

    def x25519(self, k_bytes: bytes, u_bytes: bytes) -> bytes:
      """
      Perform X25519 scalar multiplication using the optimized Montgomery ladder.

      Steps:
        1. Clamp the 32-byte scalar (per RFC 7748).
        2. Decode the base points x-coordinate from u_bytes.
        3. Run the optimized ladder to compute the scalar multiple.
        4. Encode the resulting x-coordinate as 32 bytes.
      """
      k_int = clamp_scalar(bytearray(k_bytes))
      xP = decode_u(u_bytes)
      result_int = self.ladder(k_int, xP)
      return encode_u_coordinate(result_int)