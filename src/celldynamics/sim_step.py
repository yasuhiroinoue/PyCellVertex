from __future__ import annotations

from .constants import DELTA_TIME
from .force import calc_area_force, calc_line_force, reset_forces, reset_line_tension
from .model import GlobalState


def update_pulse(state: GlobalState, delta_time: float = DELTA_TIME) -> None:
    for cp in state.p_c:
        cp.cell_time += delta_time


def compute_midpoint_stage1_forces(state: GlobalState, power_pcp: float = 2.0) -> None:
    reset_forces(state, deg=0)
    reset_line_tension(state)
    calc_line_force(state, deg=0, power_pcp=power_pcp)
    calc_area_force(state, deg=0)


def motion_vertex_second_step(
    state: GlobalState,
    power_pcp: float = 2.0,
    delta_time: float = DELTA_TIME,
) -> None:
    """Python port of C++ ODE_solver::motionVertexSecond for one step."""
    # stage1 forces on loc[0]
    reset_forces(state, deg=0)
    reset_line_tension(state)
    calc_line_force(state, deg=0, power_pcp=power_pcp)
    calc_area_force(state, deg=0)

    for vp in state.p_v:
        vp.loc[1] = vp.loc[0] + vp.frc[0] * (delta_time / 2.0)
        vp.frc[0] = vp.frc[0] * 0.0

    # stage2 forces on loc[1]
    reset_forces(state, deg=1)
    calc_line_force(state, deg=1, power_pcp=power_pcp)
    calc_area_force(state, deg=1)

    for vp in state.p_v:
        vp.loc[0] += vp.frc[1] * delta_time
        vp.frc[1] = vp.frc[1] * 0.0
