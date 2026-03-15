from __future__ import annotations

import math
import numpy as np

from .geometry import area_gradient, polygon_signed_area
from .model import GlobalState
from .vec import Vec3


def calc_area_force(state: GlobalState, deg: int) -> None:
    """Vectorized calculation of area elasticity force using NumPy."""
    num_v = len(state.p_v)
    num_c = len(state.p_c)
    
    if num_c == 0:
        return
        
    # --- 1. Extract data ---
    pos = np.array([[vp.loc[deg].x, vp.loc[deg].y, vp.loc[deg].z] for vp in state.p_v], dtype=np.float64)
    
    counts = np.array([len(cp.vi) for cp in state.p_c], dtype=np.int32)
    starts = np.zeros(num_c, dtype=np.int32)
    if num_c > 1:
        starts[1:] = np.cumsum(counts)[:-1]
    
    flat_vi = np.array([vid for cp in state.p_c for vid in cp.vi], dtype=np.int32)
    if len(flat_vi) == 0:
        return
        
    k_area = np.array([cp.K_AREA for cp in state.p_c], dtype=np.float64)
    area_eq = np.array([cp.AREA_EQ for cp in state.p_c], dtype=np.float64)
    
    # --- 2. Calculate next/prev indices for polygon edges ---
    total_vi = len(flat_vi)
    counts_rep = np.repeat(counts, counts)
    starts_rep = np.repeat(starts, counts)
    local_idx = np.arange(total_vi, dtype=np.int32) - starts_rep
    
    next_local = (local_idx + 1) % counts_rep
    prev_local = (local_idx - 1) % counts_rep
    
    flat_next = starts_rep + next_local
    flat_prev = starts_rep + prev_local
    
    # --- 3. Gather coordinates ---
    v_curr = pos[flat_vi]
    v_next = pos[flat_vi[flat_next]]
    v_prev = pos[flat_vi[flat_prev]]
    
    # --- 4. Calculate areas ---
    cross_z = 0.5 * (v_curr[:, 0] * v_next[:, 1] - v_next[:, 0] * v_curr[:, 1])
    
    cell_indices = np.repeat(np.arange(num_c, dtype=np.int32), counts)
    cell_area = np.zeros(num_c, dtype=np.float64)
    np.add.at(cell_area, cell_indices, cross_z)
    
    factor = -1.0 * k_area * (cell_area - area_eq)
    factor_rep = factor[cell_indices]
    
    # --- 5. Calculate area gradients and forces ---
    grad_x = 0.5 * (v_next[:, 1] - v_prev[:, 1])
    grad_y = 0.5 * (v_prev[:, 0] - v_next[:, 0])
    
    frc_x = factor_rep * grad_x
    frc_y = factor_rep * grad_y
    
    # --- 6. Aggregate forces to vertices and write back ---
    vertex_frc_x = np.zeros(num_v, dtype=np.float64)
    vertex_frc_y = np.zeros(num_v, dtype=np.float64)
    
    np.add.at(vertex_frc_x, flat_vi, frc_x)
    np.add.at(vertex_frc_y, flat_vi, frc_y)
    
    for i, vp in enumerate(state.p_v):
        vp.frc[deg].x += float(vertex_frc_x[i])
        vp.frc[deg].y += float(vertex_frc_y[i])


def calc_line_force(state: GlobalState, deg: int, power_pcp: float, delta_time: float = 1e-4) -> None:
    """Vectorized calculation of line tension and PCP forces using NumPy."""
    num_v = len(state.p_v)
    num_l = len(state.p_l)
    num_c = len(state.p_c)
    
    if num_l == 0:
        return

    # --- 1. Extract data ---
    pos = np.array([[vp.loc[deg].x, vp.loc[deg].y, vp.loc[deg].z] for vp in state.p_v], dtype=np.float64)
    l_vi = np.array([lp.vi for lp in state.p_l], dtype=np.int32)
    
    k2 = np.array([lp.K2_LENGTH for lp in state.p_l], dtype=np.float64)
    k1 = np.array([lp.K1_LENGTH for lp in state.p_l], dtype=np.float64)
    k1_pcp = np.array([lp.K1_PCP_LENGTH for lp in state.p_l], dtype=np.float64)
    length_eq = np.array([lp.LENGTH_EQ for lp in state.p_l], dtype=np.float64)
    
    # --- 2. Edge vectors and lengths ---
    v0_idx = l_vi[:, 0]
    v1_idx = l_vi[:, 1]
    
    edge_vec = pos[v0_idx] - pos[v1_idx]
    edge_len = np.linalg.norm(edge_vec, axis=1)
    
    l_grad = np.zeros_like(edge_vec)
    valid_len = edge_len > 1e-5
    l_grad[valid_len] = edge_vec[valid_len] / edge_len[valid_len, np.newaxis]
    
    # --- 3. Cell-to-Line mappings ---
    line_to_cell_lidx = []
    line_to_cell_cidx = []
    for lidx, lp in enumerate(state.p_l):
        for cidx in lp.ci:
            if cidx != -1:
                line_to_cell_lidx.append(lidx)
                line_to_cell_cidx.append(cidx)
                
    map_lidx = np.array(line_to_cell_lidx, dtype=np.int32)
    map_cidx = np.array(line_to_cell_cidx, dtype=np.int32)
    
    # --- 4. K2 force calculation ---
    cell_i_length = np.zeros(num_c, dtype=np.float64)
    np.add.at(cell_i_length, map_cidx, edge_len[map_lidx])
    
    i_len_diff = cell_i_length[map_cidx] - length_eq[map_lidx]
    edge_k2_factor = np.zeros(num_l, dtype=np.float64)
    np.add.at(edge_k2_factor, map_lidx, i_len_diff)
    edge_k2_factor *= k2
    
    frc_k2 = l_grad * edge_k2_factor[:, np.newaxis]
    
    # --- 5. K1 and PCP force calculation ---
    div = l_grad[:, 0]
    
    t_offset = (delta_time * 0.5) if deg == 1 else 0.0
    cell_time = np.array([cp.cell_time for cp in state.p_c], dtype=np.float64)
    cell_T = np.array([cp.cell_T for cp in state.p_c], dtype=np.float64)
    cell_phase = np.array([cp.cell_phase for cp in state.p_c], dtype=np.float64)
    
    t_eval = cell_time + t_offset
    sint_t_cell = np.sin(2.0 * np.pi * t_eval / cell_T - cell_phase)
    
    mapped_div = div[map_lidx]
    mapped_k1 = k1[map_lidx]
    mapped_k1_pcp = k1_pcp[map_lidx]
    mapped_sint_t = sint_t_cell[map_cidx]
    
    # Use np.abs for negative div when power is float, or since power is typically integer (2.0), np.power is fine.
    # original code: div**power_pcp
    # Using np.sign and np.abs can handle weird negative fraction powers if they exist, but normally it's 2.
    mapped_pcp_tension = (np.abs(mapped_div) ** power_pcp) * np.sign(mapped_div)**int(power_pcp) * mapped_k1_pcp * (mapped_sint_t + 1.0)
    
    edge_k1_pcp_factor = np.zeros(num_l, dtype=np.float64)
    np.add.at(edge_k1_pcp_factor, map_lidx, mapped_k1 + mapped_pcp_tension)
    
    edge_lt_update = np.zeros(num_l, dtype=np.float64)
    np.add.at(edge_lt_update, map_lidx, mapped_pcp_tension)
    
    frc_k1_pcp = l_grad * edge_k1_pcp_factor[:, np.newaxis]
    
    # --- 6. Aggregate forces to vertices ---
    total_edge_frc = frc_k2 + frc_k1_pcp
    vertex_frc = np.zeros((num_v, 3), dtype=np.float64)
    
    np.add.at(vertex_frc, v0_idx, -total_edge_frc)
    np.add.at(vertex_frc, v1_idx, total_edge_frc)
    
    # --- 7. Write back to Python objects ---
    for i, vp in enumerate(state.p_v):
        vp.frc[deg].x += float(vertex_frc[i, 0])
        vp.frc[deg].y += float(vertex_frc[i, 1])
        vp.frc[deg].z += float(vertex_frc[i, 2])
        
    for i, lp in enumerate(state.p_l):
        lp.lt += float(edge_lt_update[i])


def omp_reduction_frc(state: GlobalState, deg: int) -> None:
    for vp in state.p_v:
        vp.frc[deg] = Vec3(0.0, 0.0, 0.0)
        for tnum in range(state.thread_num):
            vp.frc[deg] += vp.frc_thread[tnum]
            vp.frc_thread[tnum] = Vec3(0.0, 0.0, 0.0)


def omp_reduction_lt(state: GlobalState) -> None:
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
