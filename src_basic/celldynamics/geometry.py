from __future__ import annotations

from collections.abc import Sequence

from .vec import Vec3


def polygon_signed_area(vertices: Sequence[Vec3]) -> float:
    """Signed area in xy plane, parity with C++ force::calcAreaForce."""
    n = len(vertices)
    area = 0.0
    for i in range(n):
        v0 = vertices[i]
        v1 = vertices[(i + 1) % n]
        area += 0.5 * (v0.x * v1.y - v1.x * v0.y)
    return area


def polygon_perimeter(vertices: Sequence[Vec3]) -> float:
    n = len(vertices)
    length = 0.0
    for i in range(n):
        edge = vertices[(i + 1) % n] - vertices[i]
        length += edge.norm()
    return length


def area_gradient(vertices: Sequence[Vec3], index: int) -> Vec3:
    """Area gradient for one vertex, equivalent to C++ `s_grad` logic."""
    n = len(vertices)
    prev_v = vertices[(index - 1) % n]
    next_v = vertices[(index + 1) % n]
    return Vec3(
        0.5 * (next_v.y - prev_v.y),
        0.5 * (prev_v.x - next_v.x),
        0.0,
    )
