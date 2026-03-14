from __future__ import annotations

from dataclasses import dataclass, field

from .vec import Vec3


@dataclass(slots=True)
class Vertex:
    """Minimal vertex structure for force module parity tests."""

    loc: list[Vec3]
    li: list[int] = field(default_factory=list)
    ci: list[int] = field(default_factory=list)
    frc: list[Vec3] = field(default_factory=list)
    frc_thread: list[Vec3] = field(default_factory=list)

    def __post_init__(self) -> None:
        if not self.frc:
            self.frc = [Vec3() for _ in self.loc]


@dataclass(slots=True)
class Line:
    """Minimal line structure for line-force parity tests."""

    vi: tuple[int, int]
    ci: list[int]
    K1_LENGTH: float
    K1_PCP_LENGTH: float
    K2_LENGTH: float
    LENGTH_EQ: float
    lt: float = 0.0
    lt_thread: list[float] = field(default_factory=list)


@dataclass(slots=True)
class Cellula:
    """Minimal cell structure carrying topology and parameters."""

    vi: list[int]
    li: list[int]
    K_AREA: float
    AREA_EQ: float
    center: Vec3 = field(default_factory=Vec3)
    cell_time: float = 0.0
    cell_phase: float = 0.0
    cell_T: float = 1.0
    fix: int = 0


@dataclass(slots=True)
class GlobalState:
    """Minimal global container equivalent to C++ Global pointer use."""

    p_v: list[Vertex]
    p_c: list[Cellula]
    p_l: list[Line] = field(default_factory=list)
    thread_num: int = 1

    def __post_init__(self) -> None:
        for vp in self.p_v:
            if len(vp.frc_thread) != self.thread_num:
                vp.frc_thread = [Vec3() for _ in range(self.thread_num)]
        for lp in self.p_l:
            if len(lp.lt_thread) != self.thread_num:
                lp.lt_thread = [0.0 for _ in range(self.thread_num)]
