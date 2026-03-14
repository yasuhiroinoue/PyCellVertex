from __future__ import annotations

from dataclasses import replace

from .model import Cellula, GlobalState, Line, Vertex
from .vec import Vec3

EPS = 1e-9


def polygon_center(state: GlobalState, cidx: int) -> Vec3:
    cp = state.p_c[cidx]
    tmp = Vec3(0.0, 0.0, 0.0)
    area_all = 0.0
    n = len(cp.vi)
    for i in range(n):
        r1 = state.p_v[cp.vi[i]].loc[0]
        r2 = state.p_v[cp.vi[(i + 1) % n]].loc[0]
        area_tmp = 0.5 * (r1.cross(r2)).z
        area_all += area_tmp
        tmp += (r1 + r2) / 3.0 * area_tmp
    return tmp / area_all


def update_centers(state: GlobalState) -> None:
    for i, _ in enumerate(state.p_c):
        state.p_c[i].center = polygon_center(state, i)


def ccw(a: Vec3, b: Vec3, c: Vec3) -> int:
    bb = b - a
    cc = c - a
    z = bb.cross(cc).z
    if z > 0:
        return 1
    if z < 0:
        return -1
    if bb * cc < 0:
        return 2
    if bb.norm() < cc.norm():
        return -2
    return 0


def is_convex(state: GlobalState, cidx: int) -> bool:
    vi = state.p_c[cidx].vi
    n = len(vi)
    for i in range(n):
        prev = state.p_v[vi[(i - 1 + n) % n]].loc[0]
        cur = state.p_v[vi[i]].loc[0]
        nxt = state.p_v[vi[(i + 1) % n]].loc[0]
        if ccw(prev, cur, nxt) < 0:
            return False
    return True


def _intersect_segment_line(a: Vec3, b: Vec3, l0: Vec3, l1: Vec3) -> bool:
    if abs((b - a).cross(l1 - l0).z) < EPS:
        return False
    m1 = (l1 - l0).cross(a - l0).z
    m2 = (l1 - l0).cross(b - l0).z
    return m1 * m2 <= EPS


def _cross_point(a: Vec3, b: Vec3, l0: Vec3, l1: Vec3) -> Vec3:
    sv = b - a
    tv = l1 - l0
    length = tv.cross(l0 - a).z / tv.cross(sv).z
    return a + sv * length


def _erase_one(xs: list[int], x: int) -> None:
    if x in xs:
        xs.remove(x)


def sort_counter_clockwise(state: GlobalState) -> None:
    for cidx, cp in enumerate(state.p_c):
        tmp_vi: list[int] = []
        tmp_li: list[int] = []
        if not cp.vi:
            continue
        prev = -1
        cur = cp.vi[0]
        tmp_vi.append(cur)
        while len(tmp_li) < len(cp.li):
            found = False
            for lidx in cp.li:
                l = state.p_l[lidx]
                if l.vi[0] == cur and l.vi[1] != prev:
                    prev = cur
                    cur = l.vi[1]
                    tmp_li.append(lidx)
                    if cur != tmp_vi[0]:
                        tmp_vi.append(cur)
                    found = True
                    break
                if l.vi[1] == cur and l.vi[0] != prev:
                    prev = cur
                    cur = l.vi[0]
                    tmp_li.append(lidx)
                    if cur != tmp_vi[0]:
                        tmp_vi.append(cur)
                    found = True
                    break
            if not found:
                import logging
                logging.error(f"Topology breakdown in sort_counter_clockwise for cell {cidx} at vertex {cur}. Cell edges: {cp.li}, sorted so far: {tmp_li}")
                break
        area = 0.0
        if len(tmp_vi) >= 3:
            for i in range(len(tmp_vi)):
                r1 = state.p_v[tmp_vi[i]].loc[0]
                r2 = state.p_v[tmp_vi[(i + 1) % len(tmp_vi)]].loc[0]
                area += 0.5 * r1.cross(r2).z
            if area < 0:
                tmp_vi.reverse()
                tmp_li.reverse()
        cp.vi = tmp_vi
        cp.li = tmp_li


def cell_division(state: GlobalState, cellula_idx: int, axis: Vec3) -> None:
    cp = state.p_c[cellula_idx]
    center = cp.center
    l0, l1 = center, center + axis

    # anticlockwise line order for the cell
    anti_l: list[int] = []
    n = len(cp.vi)
    for i in range(n):
        a = cp.vi[i]
        b = cp.vi[(i + 1) % n]
        for lidx in cp.li:
            lv = state.p_l[lidx].vi
            if lv[0] == a and lv[1] == b:
                anti_l.append(lidx)
                break
            if lv[0] == b and lv[1] == a:
                state.p_l[lidx].vi = (a, b)
                anti_l.append(lidx)
                break

    cross: list[tuple[int, Vec3]] = []
    epoint: list[int] = []
    for lidx in anti_l:
        v0, v1 = state.p_l[lidx].vi
        a = state.p_v[v0].loc[0]
        b = state.p_v[v1].loc[0]
        if _intersect_segment_line(a, b, l0, l1):
            epoint.extend([v0, v1])
            cross.append((lidx, _cross_point(a, b, l0, l1)))

    if len(cross) != 2:
        raise RuntimeError("Cannot divide cellula")

    lim1 = lim2 = -1
    for i, vidx in enumerate(cp.vi):
        if vidx == epoint[1]:
            lim1 = i
        if vidx == epoint[2]:
            lim2 = i - 1

    v1_idx = len(state.p_v)
    v2_idx = len(state.p_v) + 1
    l1_idx = len(state.p_l)
    l2_idx = len(state.p_l) + 1
    l3_idx = len(state.p_l) + 2
    c1_idx = len(state.p_c)

    cids0 = list(state.p_l[cross[0][0]].ci)
    cids1 = list(state.p_l[cross[1][0]].ci)

    v1 = Vertex(
        loc=[cross[0][1].copy(), cross[0][1].copy()],
        li=[cross[0][0], l1_idx, l3_idx],
        ci=cids0 + [c1_idx],
        frc_thread=[Vec3() for _ in range(state.thread_num)],
    )
    v2 = Vertex(
        loc=[cross[1][1].copy(), cross[1][1].copy()],
        li=[cross[1][0], l2_idx, l3_idx],
        ci=cids1 + [c1_idx],
        frc_thread=[Vec3() for _ in range(state.thread_num)],
    )

    state.p_v[epoint[1]].li.append(l1_idx)
    _erase_one(state.p_v[epoint[1]].li, cross[0][0])
    state.p_v[epoint[1]].ci.append(c1_idx)
    _erase_one(state.p_v[epoint[1]].ci, cellula_idx)

    state.p_v[epoint[2]].li.append(l2_idx)
    _erase_one(state.p_v[epoint[2]].li, cross[1][0])
    state.p_v[epoint[2]].ci.append(c1_idx)
    _erase_one(state.p_v[epoint[2]].ci, cellula_idx)

    for i in range(lim1 + 1, lim2 + 1):
        vidx = cp.vi[i]
        state.p_v[vidx].ci.append(c1_idx)
        _erase_one(state.p_v[vidx].ci, cellula_idx)

    p0 = state.p_l[cross[0][0]]
    p1 = state.p_l[cross[1][0]]
    l1n = Line((v1_idx, epoint[1]), [x for x in p0.ci if x != cellula_idx] + [c1_idx], p0.K1_LENGTH, p0.K1_PCP_LENGTH, p0.K2_LENGTH, p0.LENGTH_EQ, lt_thread=[0.0 for _ in range(state.thread_num)])
    l2n = Line((v2_idx, epoint[2]), [x for x in p1.ci if x != cellula_idx] + [c1_idx], p1.K1_LENGTH, p1.K1_PCP_LENGTH, p1.K2_LENGTH, p1.LENGTH_EQ, lt_thread=[0.0 for _ in range(state.thread_num)])
    l3n = Line((v1_idx, v2_idx), [cellula_idx, c1_idx], p0.K1_LENGTH, p0.K1_PCP_LENGTH, p0.K2_LENGTH, p0.LENGTH_EQ, lt_thread=[0.0 for _ in range(state.thread_num)])

    state.p_l[cross[0][0]].vi = (state.p_l[cross[0][0]].vi[0], v1_idx)
    state.p_l[cross[1][0]].vi = (v2_idx, state.p_l[cross[1][0]].vi[1])

    for i in range(lim1, lim2 + 1):
        lidx = anti_l[i]
        state.p_l[lidx].ci.append(c1_idx)
        _erase_one(state.p_l[lidx].ci, cellula_idx)

    c1_vi = [v1_idx] + cp.vi[lim1 : lim2 + 2] + [v2_idx]
    c1_li = [l1_idx] + anti_l[lim1 : lim2 + 1] + [l2_idx, l3_idx]
    c1 = Cellula(vi=c1_vi, li=c1_li, K_AREA=cp.K_AREA, AREA_EQ=cp.AREA_EQ, center=cp.center.copy(), cell_time=0.0, cell_phase=0.0, cell_T=cp.cell_T, fix=0)

    cp.vi = cp.vi[:lim1] + cp.vi[lim2 + 2 :] + [v1_idx, v2_idx]
    for i in range(lim1, lim2 + 1):
        _erase_one(cp.li, anti_l[i])
    cp.li.append(l3_idx)

    # neighboring cells attached to split edges
    adj1 = next((c for c in state.p_l[cross[0][0]].ci if c != cellula_idx), -1)
    if adj1 != -1:
        state.p_c[adj1].vi.append(v1_idx)
        state.p_c[adj1].li.append(l1_idx)
    adj2 = next((c for c in state.p_l[cross[1][0]].ci if c != cellula_idx), -1)
    if adj2 != -1:
        state.p_c[adj2].vi.append(v2_idx)
        state.p_c[adj2].li.append(l2_idx)

    state.p_c.append(c1)
    state.p_l.extend([l1n, l2n, l3n])
    state.p_v.extend([v1, v2])

    sort_counter_clockwise(state)
    update_centers(state)