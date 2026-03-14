from __future__ import annotations

import math
from dataclasses import dataclass

from .division import update_centers
from .model import Cellula, GlobalState, Line, Vertex
from .vec import Vec3


@dataclass(frozen=True)
class InitParams:
    num_x: int = 5
    num_y: int = 5
    degree_accuracy: int = 2
    k_area: float = 1.0
    area_eq: float = 2.60
    k2_length: float = 1.0
    k1_length: float = -6.0
    k1_pcp_length: float = 0.0
    length_eq: float = 0.0
    pulse_t: float = 55.0
    phase_x: float = 0.7854
    phase_y: float = 1.5708


def _loc_slots(v: Vec3, degree_accuracy: int) -> list[Vec3]:
    return [v.copy() for _ in range(degree_accuracy)]


def init_plain(params: InitParams = InitParams()) -> GlobalState:
    p_v: list[Vertex] = []
    p_l: list[Line] = []
    p_c: list[Cellula] = []

    y_tmp = 0.0
    for j in range(2 * params.num_y + 1):
        for i in range(2 * params.num_x):
            if j % 2 == 0:
                if i == 0:
                    x = 0.5
                elif i % 2 == 1:
                    x = p_v[-1].loc[0].x + 1.0
                else:
                    x = p_v[-1].loc[0].x + 2.0
            else:
                if i == 0:
                    x = 0.0
                elif i % 2 == 1:
                    x = p_v[-1].loc[0].x + 2.0
                else:
                    x = p_v[-1].loc[0].x + 1.0

            p_v.append(Vertex(loc=_loc_slots(Vec3(x, y_tmp, 0.0), params.degree_accuracy)))
        y_tmp += math.sqrt(3.0) / 2.0

    for i in range(params.num_x):
        for j in range(2 * params.num_y):
            p_l.append(Line(vi=(2 * i + 2 * params.num_x * j, 2 * i + 2 * params.num_x * (j + 1)), ci=[], K1_LENGTH=params.k1_length, K1_PCP_LENGTH=params.k1_pcp_length, K2_LENGTH=params.k2_length, LENGTH_EQ=params.length_eq))
        for j in range(params.num_y + 1):
            p_l.append(Line(vi=(2 * i + 4 * params.num_x * j, 1 + 2 * i + 4 * params.num_x * j), ci=[], K1_LENGTH=params.k1_length, K1_PCP_LENGTH=params.k1_pcp_length, K2_LENGTH=params.k2_length, LENGTH_EQ=params.length_eq))
        for j in range(2 * params.num_y):
            p_l.append(Line(vi=(1 + 2 * i + 2 * params.num_x * j, 1 + 2 * i + 2 * params.num_x * (j + 1)), ci=[], K1_LENGTH=params.k1_length, K1_PCP_LENGTH=params.k1_pcp_length, K2_LENGTH=params.k2_length, LENGTH_EQ=params.length_eq))
        if i == params.num_x - 1:
            break
        for j in range(params.num_y):
            p_l.append(Line(vi=(2 * (params.num_x + 1) + 2 * i + 4 * params.num_x * j - 1, 2 * (params.num_x + 1) + 2 * i + 4 * params.num_x * j), ci=[], K1_LENGTH=params.k1_length, K1_PCP_LENGTH=params.k1_pcp_length, K2_LENGTH=params.k2_length, LENGTH_EQ=params.length_eq))

    for i in range(2 * params.num_x - 1):
        if i % 2 == 0:
            for j in range(params.num_y):
                li = [
                    2 * j + (6 * params.num_y + 1) * i // 2,
                    2 * j + 1 + (6 * params.num_y + 1) * i // 2,
                    2 * params.num_y + j + (6 * params.num_y + 1) * i // 2,
                    2 * params.num_y + (j + 1) + (6 * params.num_y + 1) * i // 2,
                    (3 * params.num_y + 1) + 2 * j + (6 * params.num_y + 1) * i // 2,
                    (3 * params.num_y + 1) + 2 * j + 1 + (6 * params.num_y + 1) * i // 2,
                ]
                vi = [
                    4 * params.num_x * j + i,
                    4 * params.num_x * j + i + 1,
                    2 * params.num_x + 4 * params.num_x * j + i + 1,
                    4 * params.num_x * (j + 1) + i + 1,
                    4 * params.num_x * (j + 1) + i,
                    2 * params.num_x + 4 * params.num_x * j + i,
                ]
                p_c.append(Cellula(vi=vi, li=li, K_AREA=params.k_area, AREA_EQ=params.area_eq))
        else:
            for j in range(params.num_y - 1):
                li = [
                    2 * j + (3 * params.num_y + 1) + 1 + (6 * params.num_y + 1) * (i - 1) // 2,
                    2 * j + (3 * params.num_y + 1) + 2 + (6 * params.num_y + 1) * (i - 1) // 2,
                    j + (5 * params.num_y + 1) + (6 * params.num_y + 1) * (i - 1) // 2,
                    (j + 1) + (5 * params.num_y + 1) + (6 * params.num_y + 1) * (i - 1) // 2,
                    2 * j + 1 + (6 * params.num_y + 1) * (i + 1) // 2,
                    2 * j + 2 + (6 * params.num_y + 1) * (i + 1) // 2,
                ]
                vi = [
                    2 * params.num_x + 1 + 4 * params.num_x * j + (i - 1),
                    2 * params.num_x + 2 + 4 * params.num_x * j + (i - 1),
                    4 * params.num_x + 2 + 4 * params.num_x * j + (i - 1),
                    6 * params.num_x + 2 + 4 * params.num_x * j + (i - 1),
                    6 * params.num_x + 1 + 4 * params.num_x * j + (i - 1),
                    4 * params.num_x + 1 + 4 * params.num_x * j + (i - 1),
                ]
                p_c.append(Cellula(vi=vi, li=li, K_AREA=params.k_area, AREA_EQ=params.area_eq))

    # vertex-line incidence
    for lidx, lp in enumerate(p_l):
        p_v[lp.vi[0]].li.append(lidx)
        p_v[lp.vi[1]].li.append(lidx)

    # line-cell and vertex-cell incidence
    for cidx, cp in enumerate(p_c):
        for lidx in cp.li:
            p_l[lidx].ci.append(cidx)
        for vidx in cp.vi:
            p_v[vidx].ci.append(cidx)

    flag = 0
    flag_l = 0
    flag_s = 0
    interval = params.num_y
    phase_tmp = 0.0
    phase_x_tmp = 0.0

    for i, cp in enumerate(p_c):
        count = (i + 1) - flag_l * params.num_y - flag_s * (params.num_y - 1)
        if count % interval == 0:
            if flag == 0:
                flag = 1
                flag_l += 1
                interval = params.num_y - 1
                phase_x_tmp += params.phase_x
                phase_tmp = phase_x_tmp
            else:
                flag = 0
                flag_s += 1
                interval = params.num_y
                phase_x_tmp += params.phase_x
                phase_x_tmp = phase_x_tmp
        cp.cell_time = 0.0
        cp.cell_phase = phase_tmp
        cp.cell_T = params.pulse_t
        phase_tmp += params.phase_y

    state = GlobalState(p_v=p_v, p_c=p_c, p_l=p_l)
    update_centers(state)
    return state
