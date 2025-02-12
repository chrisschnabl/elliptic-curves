from abc import ABC, abstractmethod
from util import clamp_scalar, decode_u, modinv, encode_u_coordinate
from curve import AffinePoint, Curve, IdentityPoint


class X25519Curve(Curve):
    """
    An X25519 is a curve of form X25519, that provides the x25519 method
    It takes bytes and returns bytes accoridng to the RFC7748 and uses some scalar_mult method
    """
    def __init__(self):
       self.p = 2**255 - 19
       self.a24 = 121665
 
    def decode_public_key(self, public_key_bytes: bytes) -> int:
        return decode_u(public_key_bytes)
    
    def recover_point(self, x: int) -> AffinePoint:
        return AffinePoint(x, 0)
    
    # For fun we could also transform between X25519 and Edwards25519
    def x25519(self, private_key_bytes: bytes, public_key_bytes: bytes) -> bytes:
      """
      Perform X25519 scalar multiplication using the optimized Montgomery ladder.

      Steps:
        1. Clamp the 32-byte scalar (per RFC 7748).
        2. Decode the base points x-coordinate from u_bytes.
        3. Run the optimized ladder to compute the scalar multiple.
        4. Encode the resulting x-coordinate as 32 bytes.
      """
      k_int = clamp_scalar(bytearray(private_key_bytes))
      xP = self.decode_public_key(public_key_bytes)
      point = self.recover_point(xP)  # Will be (Xp, 0) for Montgomery ladders
      result = self.scalar_mult(point, k_int)
      if result is IdentityPoint:
         raise ValueError("Resulting point is the point at infinity")
      return encode_u_coordinate(result.x)