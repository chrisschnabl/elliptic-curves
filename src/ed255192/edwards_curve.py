from dataclasses import dataclass
from typing import TypeVar

from util import modinv
from abc import ABC, abstractmethod

TPoint = TypeVar('TPoint', bound='AffinePoint')


@dataclass
class AffinePoint():
    x: int
    y: int

@dataclass
class ExtendedPoint(AffinePoint):
    z: int
    t: int

IdentityPoint = None
Point = AffinePoint | IdentityPoint

class EdwardsCurve(ABC):
    """
    Abstract interface for a point on an Edwards curve.
    Every (non-identity) point supports addition, scalar multiplication, and doubling.
    Identity is denoted by None.
    """
    def __init__(self) -> None:
        self.a = -1
        self.p = 2**255 - 19
        # The Edwards curve constant for Ed25519.
        self.d = (-121665 * modinv(121666, self.p)) % self.p
        # The subgroup order for Ed25519.
        self.q = 2**252 + 27742317777372353535851937790883648493
        self.B = AffinePoint(
            15112221349535400772501151409588531511454012693041857206046113283949847762202,
            46316835694926478169428394003475163141307993866256225615783033603165251855960,
        )


    @abstractmethod
    def add(self, Q: Point, R: Point) -> Point:
        raise NotImplementedError

    def scalar_mult(self, R: Point, scalar: int) -> Point:
        if R is IdentityPoint:
            return IdentityPoint

        Q = None
        while scalar > 0:
            if scalar % 2 == 1:
                Q = R if Q is None else self.add(Q, R)
            R = self.double(R)
            scalar //= 2
        return Q

    @abstractmethod
    def double(self, R: Point) -> Point:
        # Terrible default implementation
        # Often benefits from optimization 
        return self.add(R, R)
    
    @abstractmethod
    def compress(self, point: Point) -> bytes:
        raise NotImplementedError

    @abstractmethod
    def uncompress(self, comp: bytes) -> Point:
        raise NotImplementedError
    
    @abstractmethod
    def point_equals(self, P: Point, Q: Point) -> bool:
        raise NotImplementedError