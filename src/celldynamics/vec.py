from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(slots=True)
class Vec3:
    """C++ `vec.h` parity-oriented 3D vector utility."""

    x: float = 0.0
    y: float = 0.0
    z: float = 0.0

    def copy(self) -> Vec3:
        return Vec3(self.x, self.y, self.z)

    def in_(self, x: float, y: float, z: float) -> None:
        self.x = x
        self.y = y
        self.z = z

    def in_vec(self, v: Vec3) -> None:
        self.x = v.x
        self.y = v.y
        self.z = v.z

    def icast(self) -> Vec3:
        # C++ cast to int truncates toward zero.
        return Vec3(int(self.x), int(self.y), int(self.z))

    def norm(self) -> float:
        return math.sqrt(self.sqr())

    def sqr(self) -> float:
        return self.x * self.x + self.y * self.y + self.z * self.z

    def dot(self, other: Vec3) -> float:
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other: Vec3) -> Vec3:
        return Vec3(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )

    def __add__(self, other: Vec3 | float) -> Vec3:
        if isinstance(other, Vec3):
            return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)
        return Vec3(self.x + other, self.y + other, self.z + other)

    def __radd__(self, other: float) -> Vec3:
        return self + other

    def __sub__(self, other: Vec3 | float) -> Vec3:
        if isinstance(other, Vec3):
            return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)
        return Vec3(self.x - other, self.y - other, self.z - other)

    def __rsub__(self, other: float) -> Vec3:
        return Vec3(other - self.x, other - self.y, other - self.z)

    def __mul__(self, other: Vec3 | float) -> Vec3 | float:
        if isinstance(other, Vec3):
            return self.dot(other)
        return Vec3(self.x * other, self.y * other, self.z * other)

    def __rmul__(self, other: float) -> Vec3:
        return self * other

    def __truediv__(self, other: Vec3 | float) -> Vec3:
        if isinstance(other, Vec3):
            return Vec3(self.x / other.x, self.y / other.y, self.z / other.z)
        return Vec3(self.x / other, self.y / other, self.z / other)

    def __rtruediv__(self, other: float) -> Vec3:
        return Vec3(other / self.x, other / self.y, other / self.z)

    def __iadd__(self, other: Vec3) -> Vec3:
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self

    def __isub__(self, other: Vec3) -> Vec3:
        self.x -= other.x
        self.y -= other.y
        self.z -= other.z
        return self

    def __imul__(self, other: Vec3 | float) -> Vec3:
        if isinstance(other, Vec3):
            self.x *= other.x
            self.y *= other.y
            self.z *= other.z
        else:
            self.x *= other
            self.y *= other
            self.z *= other
        return self

    def __itruediv__(self, other: Vec3 | float) -> Vec3:
        if isinstance(other, Vec3):
            self.x /= other.x
            self.y /= other.y
            self.z /= other.z
        else:
            self.x /= other
            self.y /= other
            self.z /= other
        return self
