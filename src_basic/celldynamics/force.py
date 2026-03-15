from __future__ import annotations

import math

from .geometry import area_gradient, polygon_signed_area
from .model import GlobalState
from .vec import Vec3


def calc_area_force(state: GlobalState, deg: int) -> None:
    """Python port of C++ `force::calcAreaForce` (single-threaded baseline)."""
    for cp in state.p_c:
        vertices = [state.p_v[idx].loc[deg] for idx in cp.vi]
        i_area = polygon_signed_area(vertices)

        factor = -1.0 * cp.K_AREA * (i_area - cp.AREA_EQ)
        for local_idx, global_vidx in enumerate(cp.vi):
            s_grad = area_gradient(vertices, local_idx)
            frc_tmp = factor * s_grad
            state.p_v[global_vidx].frc[deg] += frc_tmp


def calc_line_force(state: GlobalState, deg: int, power_pcp: float, delta_time: float = 1e-4) -> None:
    """Python port of C++ `force::calcLineForce` (single-threaded baseline)."""
    for cp in state.p_c:
        i_length = 0.0
        for line_idx in cp.li:
            lp = state.p_l[line_idx]
            v0 = state.p_v[lp.vi[0]].loc[deg]
            v1 = state.p_v[lp.vi[1]].loc[deg]
            i_length += (v1 - v0).norm()

        for line_idx in cp.li:
            lp = state.p_l[line_idx]
            i0, i1 = lp.vi
            v0 = state.p_v[i0].loc[deg]
            v1 = state.p_v[i1].loc[deg]

            edge = v0 - v1
            edge_length = edge.norm()
            l_grad = Vec3(0.0, 0.0, 0.0)
            if edge_length > 1e-5:
                l_grad = edge / edge_length

            frc_tmp = l_grad * lp.K2_LENGTH * (i_length - lp.LENGTH_EQ)
            state.p_v[i0].frc[deg] -= frc_tmp
            state.p_v[i1].frc[deg] += frc_tmp

        for line_idx in cp.li:
            lp = state.p_l[line_idx]
            i0, i1 = lp.vi
            v0 = state.p_v[i0].loc[deg]
            v1 = state.p_v[i1].loc[deg]

            edge = v0 - v1
            edge_length = edge.norm()
            l_grad = Vec3(0.0, 0.0, 0.0)
            if edge_length > 1e-5:
                l_grad = edge / edge_length

            frc_tmp = l_grad * lp.K1_LENGTH

            pcp_vec = Vec3(1.0, 0.0, 0.0)
            pcp_vec /= pcp_vec.norm()
            div = l_grad * pcp_vec
            
            t_offset = (delta_time * 0.5) if deg == 1 else 0.0
            t_eval = cp.cell_time + t_offset
            sint_t = math.sin(2.0 * math.pi * t_eval / cp.cell_T - cp.cell_phase)
            
            line_tension_pcp = (div**power_pcp) * lp.K1_PCP_LENGTH * (sint_t + 1.0)

            frc_tmp += l_grad * line_tension_pcp
            lp.lt += line_tension_pcp

            state.p_v[i0].frc[deg] -= frc_tmp
            state.p_v[i1].frc[deg] += frc_tmp


def omp_reduction_frc(state: GlobalState, deg: int) -> None:
    """Parity helper for C++ `OMP_Reduction_Frc` without stochastic fluctuation."""
    for vp in state.p_v:
        vp.frc[deg] = Vec3(0.0, 0.0, 0.0)
        for tnum in range(state.thread_num):
            vp.frc[deg] += vp.frc_thread[tnum]
            vp.frc_thread[tnum] = Vec3(0.0, 0.0, 0.0)


def omp_reduction_lt(state: GlobalState) -> None:
    """Parity helper for C++ `OMP_Reduction_Lt`."""
    for lp in state.p_l:
        lp.lt = 0.0
        for tnum in range(state.thread_num):
            lp.lt += lp.lt_thread[tnum]
            lp.lt_thread[tnum] = 0.0


def reset_forces(state: GlobalState, deg: int) -> None:
    for vp in state.p_v:
        vp.frc[deg] = Vec3(0.0, 0.0, 0.0)
        for tnum in range(state.thread_num):
            vp.frc_thread[tnum] = Vec3(0.0, 0.0, 0.0)


def reset_line_tension(state: GlobalState) -> None:
    for lp in state.p_l:
        lp.lt = 0.0
        for tnum in range(state.thread_num):
            lp.lt_thread[tnum] = 0.0
