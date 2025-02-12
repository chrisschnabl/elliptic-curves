from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import TypeVar, override

from util import decode_u


TPoint = TypeVar('TPoint', bound='AffinePoint')


@dataclass
class AffinePoint():
    x: int
    y: int

IdentityPoint = None
Point = AffinePoint | IdentityPoint

# TODO: unify this with the x225519 curve which this is equal to
class Curve(ABC):
    @abstractmethod
    def scalar_mult(self, R: Point, scalar: int) -> Point:
        raise NotImplementedError


class DoubleAndAddCurve(Curve):
    @abstractmethod
    def add(self, Q: Point, R: Point) -> Point:
        raise NotImplementedError

    @override
    def scalar_mult(self, R: Point, scalar: int) -> Point:
        if R is IdentityPoint:
            return IdentityPoint

        Q = None
        while scalar > 0:
            if scalar % 2 == 1:
                # TODO: check if this also works for x25519
                Q = R if Q is None else self.add(Q, R)
            R = self.double(R)
            scalar //= 2
        return Q

    @abstractmethod
    def double(self, R: Point) -> Point:
        # Terrible default implementation
        # Often benefits from optimization 
        return self.add(R, R)
    