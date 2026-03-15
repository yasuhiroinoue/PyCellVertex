from __future__ import annotations

import numpy as np

from .constants import DELTA_TIME
from .force import calc_area_force, calc_line_force, reset_forces, reset_line_tension
from .model import GlobalState


def update_pulse(state: GlobalState, delta_time: float = DELTA_TIME) -> None:
    for cp in state.p_c:
        cp.cell_time += delta_time


def compute_midpoint_stage1_forces(
    state: GlobalState, power_pcp: float = 2.0, delta_time: float = DELTA_TIME
) -> None:
    reset_forces(state, deg=0)
    reset_line_tension(state)
    calc_line_force(state, deg=0, power_pcp=power_pcp, delta_time=delta_time)
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
    calc_line_force(state, deg=0, power_pcp=power_pcp, delta_time=delta_time)
    calc_area_force(state, deg=0)

    n_v = len(state.p_v)
    if n_v > 0:
        loc0 = np.empty((n_v, 3), dtype=np.float64)
        frc0 = np.empty((n_v, 3), dtype=np.float64)

        for i, vp in enumerate(state.p_v):
            loc0[i, 0] = vp.loc[0].x
            loc0[i, 1] = vp.loc[0].y
            loc0[i, 2] = vp.loc[0].z
            frc0[i, 0] = vp.frc[0].x
            frc0[i, 1] = vp.frc[0].y
            frc0[i, 2] = vp.frc[0].z

        loc1 = loc0 + frc0 * (delta_time / 2.0)

        for i, vp in enumerate(state.p_v):
            vp.loc[1].in_(loc1[i, 0], loc1[i, 1], loc1[i, 2])
            vp.frc[0].in_(0.0, 0.0, 0.0)

    # stage2 forces on loc[1]
    reset_forces(state, deg=1)
    calc_line_force(state, deg=1, power_pcp=power_pcp, delta_time=delta_time)
    calc_area_force(state, deg=1)

    if n_v > 0:
        loc0 = np.empty((n_v, 3), dtype=np.float64)
        frc1 = np.empty((n_v, 3), dtype=np.float64)

        for i, vp in enumerate(state.p_v):
            loc0[i, 0] = vp.loc[0].x
            loc0[i, 1] = vp.loc[0].y
            loc0[i, 2] = vp.loc[0].z
            frc1[i, 0] = vp.frc[1].x
            frc1[i, 1] = vp.frc[1].y
            frc1[i, 2] = vp.frc[1].z

        loc0_new = loc0 + frc1 * delta_time

        for i, vp in enumerate(state.p_v):
            vp.loc[0].in_(loc0_new[i, 0], loc0_new[i, 1], loc0_new[i, 2])
            vp.frc[1].in_(0.0, 0.0, 0.0)
