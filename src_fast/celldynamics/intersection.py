from __future__ import annotations

from .division import sort_counter_clockwise
from .model import GlobalState, Line, Vertex
from .rearrange import _erase_one
from .vec import Vec3


def _contain_vertex_in_cell(state: GlobalState, vidx: int, cidx: int) -> bool:
    ret = False
    vx = state.p_v[vidx].loc[0]
    vi = state.p_c[cidx].vi
    n = len(vi)
    for i in range(n):
        p1 = state.p_v[vi[i]].loc[0] - vx
        p2 = state.p_v[vi[(i + 1) % n]].loc[0] - vx
        if p1.y > p2.y:
            p1, p2 = p2, p1
        if p1.y <= 0 < p2.y and p1.cross(p2).z < 0:
            ret = not ret
    return ret


def _projection(a: Vec3, b: Vec3, p: Vec3) -> Vec3:
    t = ((p - a) * (a - b)) / (a - b).norm()
    return a + (a - b) * t


def _distance_line_point(a: Vec3, b: Vec3, p: Vec3) -> float:
    return (p - _projection(a, b, p)).norm()


def _closest_line_index(state: GlobalState, vidx: int, cidx: int) -> int:
    best = -1.0
    best_l = state.p_c[cidx].li[0]
    p = state.p_v[vidx].loc[0]
    for lidx in state.p_c[cidx].li:
        v0, v1 = state.p_l[lidx].vi
        d = _distance_line_point(state.p_v[v0].loc[0], state.p_v[v1].loc[0], p)
        if best < 0 or d < best:
            best = d
            best_l = lidx
    return best_l


def _cross_point(a1: Vec3, a2: Vec3, b1: Vec3, b2: Vec3) -> Vec3:
    sv = a2 - a1
    tv = b2 - b1
    length = tv.cross(b1 - a1).z / tv.cross(sv).z
    return a1 + sv * length


def _intersection_type220(
    state: GlobalState,
    vinter_idx: int,
    cobs_idx: int,
    lobs_idx: int,
    step: int | None = None,
    event_log: list[str] | None = None,
) -> None:
    vobs1_idx, vobs2_idx = state.p_l[lobs_idx].vi

    epoint: list[tuple[int, int]] = []
    for lidx in state.p_v[vinter_idx].li:
        a, b = state.p_l[lidx].vi
        if a != vinter_idx:
            epoint.append((a, lidx))
        if b != vinter_idx:
            epoint.append((b, lidx))

    vadj1, l_out1 = epoint[0]
    vadj2, l_out2 = epoint[1]

    vinter_vec = state.p_v[vinter_idx].loc[0]
    out1_vec = state.p_v[vadj1].loc[0]
    out2_vec = state.p_v[vadj2].loc[0]
    vobs1_vec = state.p_v[vobs1_idx].loc[0]
    vobs2_vec = state.p_v[vobs2_idx].loc[0]

    cp1 = _cross_point(vobs1_vec, vobs2_vec, vinter_vec, out1_vec)
    cp2 = _cross_point(vobs1_vec, vobs2_vec, vinter_vec, out2_vec)

    dist1 = (cp1 - vobs1_vec).norm()
    dist2 = (cp2 - vobs1_vec).norm()

    if dist1 > dist2:
        vadj1, vadj2 = vadj2, vadj1
        l_out1, l_out2 = l_out2, l_out1
        cp1, cp2 = cp2, cp1

    v1_idx = len(state.p_v)
    v2_idx = len(state.p_v) + 1
    l1_idx = len(state.p_l)
    l_inner1_idx = len(state.p_l) + 1
    l_inner2_idx = len(state.p_l) + 2

    cinter1_ci = state.p_l[l_out1].ci.copy()
    cinter2_ci = state.p_l[l_out2].ci.copy()

    v1_ci = list(set(cinter1_ci + [cobs_idx]))
    v2_ci = list(set(cinter2_ci + [cobs_idx]))

    v1 = Vertex(
        loc=[cp1.copy(), cp1.copy()],
        li=[l_out1, l1_idx, l_inner1_idx],
        ci=v1_ci,
        frc_thread=[Vec3() for _ in range(state.thread_num)],
    )

    v2 = Vertex(
        loc=[cp2.copy(), cp2.copy()],
        li=[l_out2, l_inner2_idx, lobs_idx],
        ci=v2_ci,
        frc_thread=[Vec3() for _ in range(state.thread_num)],
    )

    l1 = Line(
        vi=(vobs1_idx, v1_idx),
        ci=state.p_l[lobs_idx].ci.copy(),
        K1_LENGTH=state.p_l[lobs_idx].K1_LENGTH,
        K1_PCP_LENGTH=state.p_l[lobs_idx].K1_PCP_LENGTH,
        K2_LENGTH=state.p_l[lobs_idx].K2_LENGTH,
        LENGTH_EQ=state.p_l[lobs_idx].LENGTH_EQ,
        lt=0.0,
        lt_thread=[0.0 for _ in range(state.thread_num)],
    )

    l_inner1 = Line(
        vi=(v1_idx, vinter_idx),
        ci=v1_ci.copy(),
        K1_LENGTH=state.p_l[lobs_idx].K1_LENGTH,
        K1_PCP_LENGTH=state.p_l[lobs_idx].K1_PCP_LENGTH,
        K2_LENGTH=state.p_l[lobs_idx].K2_LENGTH,
        LENGTH_EQ=state.p_l[lobs_idx].LENGTH_EQ,
        lt=0.0,
        lt_thread=[0.0 for _ in range(state.thread_num)],
    )

    l_inner2 = Line(
        vi=(vinter_idx, v2_idx),
        ci=v2_ci.copy(),
        K1_LENGTH=state.p_l[lobs_idx].K1_LENGTH,
        K1_PCP_LENGTH=state.p_l[lobs_idx].K1_PCP_LENGTH,
        K2_LENGTH=state.p_l[lobs_idx].K2_LENGTH,
        LENGTH_EQ=state.p_l[lobs_idx].LENGTH_EQ,
        lt=0.0,
        lt_thread=[0.0 for _ in range(state.thread_num)],
    )

    a, b = state.p_l[lobs_idx].vi
    state.p_l[lobs_idx].vi = (v2_idx if a == vobs1_idx else a, v2_idx if b == vobs1_idx else b)

    a, b = state.p_l[l_out1].vi
    state.p_l[l_out1].vi = (v1_idx if a == vinter_idx else a, v1_idx if b == vinter_idx else b)

    a, b = state.p_l[l_out2].vi
    state.p_l[l_out2].vi = (v2_idx if a == vinter_idx else a, v2_idx if b == vinter_idx else b)

    _erase_one(state.p_v[vinter_idx].li, l_out1)
    _erase_one(state.p_v[vinter_idx].li, l_out2)
    state.p_v[vinter_idx].li.extend([l_inner1_idx, l_inner2_idx])
    if cobs_idx not in state.p_v[vinter_idx].ci:
        state.p_v[vinter_idx].ci.append(cobs_idx)

    _erase_one(state.p_v[vobs1_idx].li, lobs_idx)
    state.p_v[vobs1_idx].li.append(l1_idx)

    if vinter_idx not in state.p_c[cobs_idx].vi:
        state.p_c[cobs_idx].vi.append(vinter_idx)
    state.p_c[cobs_idx].vi.extend([v1_idx, v2_idx])
    state.p_c[cobs_idx].li.extend([l1_idx, l_inner1_idx, l_inner2_idx])

    for c in cinter1_ci:
        if v1_idx not in state.p_c[c].vi:
            state.p_c[c].vi.append(v1_idx)
        if l_inner1_idx not in state.p_c[c].li:
            state.p_c[c].li.append(l_inner1_idx)

    for c in cinter2_ci:
        if v2_idx not in state.p_c[c].vi:
            state.p_c[c].vi.append(v2_idx)
        if l_inner2_idx not in state.p_c[c].li:
            state.p_c[c].li.append(l_inner2_idx)

    state.p_v.extend([v1, v2])
    state.p_l.extend([l1, l_inner1, l_inner2])

    if event_log is not None:
        event_log.append(f"{step},intersection,220,{vinter_idx},{cobs_idx},{lobs_idx}")

def _intersection_type321(
    state: GlobalState,
    vinter_idx: int,
    cobs_idx: int,
    lobs_idx: int,
    step: int | None = None,
    event_log: list[str] | None = None,
) -> None:
    epoint: list[tuple[int, int]] = []
    for lidx in state.p_v[vinter_idx].li:
        a, b = state.p_l[lidx].vi
        if a != vinter_idx:
            epoint.append((a, lidx))
        if b != vinter_idx:
            epoint.append((b, lidx))

    vobs1 = vobs2 = linter1 = -1
    for v, l in epoint:
        a, b = state.p_l[lobs_idx].vi
        if a == v:
            vobs1, vobs2, linter1 = v, b, l
        if b == v:
            vobs1, vobs2, linter1 = v, a, l

    linter2 = vadj1 = -1
    for v, l in epoint:
        if len(state.p_l[l].ci) == 2:
            linter2, vadj1 = l, v

    linter3 = vadj2 = -1
    for v, l in epoint:
        if len(state.p_l[l].ci) == 1 and l != linter1:
            linter3, vadj2 = l, v

    if min(vobs1, vobs2, linter1, linter2, linter3, vadj1, vadj2) < 0:
        return

    cinter2 = state.p_l[linter3].ci[0]

    v1_idx = len(state.p_v)
    l1_idx = len(state.p_l)

    lobs_a, lobs_b = state.p_l[lobs_idx].vi
    cp1 = _cross_point(
        state.p_v[lobs_a].loc[0],
        state.p_v[lobs_b].loc[0],
        state.p_v[vinter_idx].loc[0],
        state.p_v[vadj2].loc[0],
    )
    cp2 = _cross_point(
        state.p_v[lobs_a].loc[0],
        state.p_v[lobs_b].loc[0],
        state.p_v[vinter_idx].loc[0],
        state.p_v[vadj1].loc[0],
    )

    vnew = Vertex(
        loc=[cp1.copy(), cp1.copy()],
        li=[lobs_idx, l1_idx, linter3],
        ci=[cobs_idx, cinter2],
        frc_thread=[Vec3() for _ in range(state.thread_num)],
    )

    _erase_one(state.p_v[vobs1].li, lobs_idx)

    _erase_one(state.p_v[vinter_idx].li, linter3)
    state.p_v[vinter_idx].li.append(l1_idx)
    state.p_v[vinter_idx].ci.append(cobs_idx)
    state.p_v[vinter_idx].loc[0] = cp2

    lnew = Line(
        vi=(v1_idx, vinter_idx),
        ci=[cobs_idx, cinter2],
        K1_LENGTH=state.p_l[lobs_idx].K1_LENGTH,
        K1_PCP_LENGTH=state.p_l[lobs_idx].K1_PCP_LENGTH,
        K2_LENGTH=state.p_l[lobs_idx].K2_LENGTH,
        LENGTH_EQ=state.p_l[lobs_idx].LENGTH_EQ,
        lt=0.0,
        lt_thread=[0.0 for _ in range(state.thread_num)],
    )

    a, b = state.p_l[lobs_idx].vi
    state.p_l[lobs_idx].vi = (v1_idx if a == vobs1 else a, v1_idx if b == vobs1 else b)
    state.p_l[linter1].ci.append(cobs_idx)

    a, b = state.p_l[linter3].vi
    state.p_l[linter3].vi = (v1_idx if a == vinter_idx else a, v1_idx if b == vinter_idx else b)

    state.p_c[cobs_idx].vi.extend([vinter_idx, v1_idx])
    state.p_c[cobs_idx].li.extend([linter1, l1_idx])

    state.p_c[cinter2].vi.append(v1_idx)
    state.p_c[cinter2].li.append(l1_idx)

    state.p_v.append(vnew)
    state.p_l.append(lnew)
    if event_log is not None:
        event_log.append(f"{step},intersection,321,{vinter_idx},{cobs_idx},{lobs_idx}")


def _intersection_type221(
    state: GlobalState,
    vinter_idx: int,
    cobs_idx: int,
    lobs_idx: int,
    step: int | None = None,
    event_log: list[str] | None = None,
) -> None:
    epoint: list[tuple[int, int]] = []
    for lidx in state.p_v[vinter_idx].li:
        a, b = state.p_l[lidx].vi
        if a != vinter_idx:
            epoint.append((a, lidx))
        if b != vinter_idx:
            epoint.append((b, lidx))

    vobs1 = vobs2 = linter1 = -1
    for v, l in epoint:
        a, b = state.p_l[lobs_idx].vi
        if a == v:
            vobs1, vobs2, linter1 = v, b, l
        if b == v:
            vobs1, vobs2, linter1 = v, a, l

    linter2 = vadj1 = -1
    for v, l in epoint:
        if v != vobs1:
            linter2, vadj1 = l, v

    if min(vobs1, vobs2, linter1, linter2, vadj1) < 0:
        return

    cp = _cross_point(
        state.p_v[vinter_idx].loc[0],
        state.p_v[vadj1].loc[0],
        state.p_v[vobs1].loc[0],
        state.p_v[vobs2].loc[0],
    )

    state.p_v[vinter_idx].li.append(lobs_idx)
    state.p_v[vinter_idx].ci.append(cobs_idx)
    state.p_v[vinter_idx].loc[0] = cp

    _erase_one(state.p_v[vobs1].li, lobs_idx)
    a, b = state.p_l[lobs_idx].vi
    state.p_l[lobs_idx].vi = (vinter_idx if a == vobs1 else a, vinter_idx if b == vobs1 else b)

    state.p_l[linter1].ci.append(cobs_idx)
    state.p_c[cobs_idx].vi.append(vinter_idx)
    state.p_c[cobs_idx].li.append(linter1)

    if event_log is not None:
        event_log.append(f"{step},intersection,221,{vinter_idx},{cobs_idx},{lobs_idx}")


def _intersection_type330(
    state: GlobalState,
    vinter_idx: int,
    cobs_idx: int,
    lobs_idx: int,
    step: int | None = None,
    event_log: list[str] | None = None,
) -> None:
    vobs1, vobs2 = state.p_l[lobs_idx].vi

    epoint: list[tuple[int, int]] = []
    for lidx in state.p_v[vinter_idx].li:
        a, b = state.p_l[lidx].vi
        if a != vinter_idx:
            epoint.append((a, lidx))
        if b != vinter_idx:
            epoint.append((b, lidx))

    linter2 = vadj2 = -1
    for v, l in epoint:
        if len(state.p_l[l].ci) == 2:
            linter2, vadj2 = l, v

    vinter = state.p_v[vinter_idx].loc[0]
    vobs1v = state.p_v[vobs1].loc[0]
    vadj2v = state.p_v[vadj2].loc[0]

    linter1 = linter3 = vadj1 = vadj3 = -1
    for v, l in epoint:
        if len(state.p_l[l].ci) == 2:
            continue
        vadjv = state.p_v[v].loc[0]
        side = (vadj2v - vinter).cross(vadjv - vinter).z * (vadj2v - vinter).cross(
            vobs1v - vinter
        ).z
        if side > 0:
            linter1, vadj1 = l, v
        if side < 0:
            linter3, vadj3 = l, v

    if min(linter1, linter2, linter3, vadj1, vadj2, vadj3) < 0:
        return

    cinter1 = state.p_l[linter1].ci[0]
    cinter2 = state.p_l[linter3].ci[0]

    cp1 = _cross_point(
        state.p_v[vobs1].loc[0], state.p_v[vobs2].loc[0], vinter, state.p_v[vadj2].loc[0]
    )
    cp2 = _cross_point(
        state.p_v[vobs1].loc[0], state.p_v[vobs2].loc[0], vinter, state.p_v[vadj3].loc[0]
    )
    cp3 = _cross_point(
        state.p_v[vobs1].loc[0], state.p_v[vobs2].loc[0], vinter, state.p_v[vadj1].loc[0]
    )

    v1_idx = len(state.p_v)
    v2_idx = len(state.p_v) + 1
    l1_idx = len(state.p_l)
    l2_idx = len(state.p_l) + 1
    l3_idx = len(state.p_l) + 2

    v1 = Vertex(
        loc=[cp1.copy(), cp1.copy()],
        li=[linter2, l2_idx, l3_idx],
        ci=[cobs_idx, cinter1, cinter2],
        frc_thread=[Vec3() for _ in range(state.thread_num)],
    )
    v2 = Vertex(
        loc=[cp2.copy(), cp2.copy()],
        li=[linter3, l3_idx, lobs_idx],
        ci=[cobs_idx, cinter2],
        frc_thread=[Vec3() for _ in range(state.thread_num)],
    )

    _erase_one(state.p_v[vinter_idx].li, linter2)
    _erase_one(state.p_v[vinter_idx].li, linter3)
    state.p_v[vinter_idx].li.extend([l1_idx, l2_idx])
    _erase_one(state.p_v[vinter_idx].ci, cinter2)
    state.p_v[vinter_idx].ci.append(cobs_idx)
    state.p_v[vinter_idx].loc[0] = cp3

    _erase_one(state.p_v[vobs1].li, lobs_idx)
    state.p_v[vobs1].li.append(l1_idx)

    l1 = Line(
        (vobs1, vinter_idx),
        [cobs_idx],
        state.p_l[lobs_idx].K1_LENGTH,
        state.p_l[lobs_idx].K1_PCP_LENGTH,
        state.p_l[lobs_idx].K2_LENGTH,
        state.p_l[lobs_idx].LENGTH_EQ,
        lt=0.0,
        lt_thread=[0.0 for _ in range(state.thread_num)],
    )
    l2 = Line(
        (vinter_idx, v1_idx),
        [cobs_idx, cinter1],
        state.p_l[lobs_idx].K1_LENGTH,
        state.p_l[lobs_idx].K1_PCP_LENGTH,
        state.p_l[lobs_idx].K2_LENGTH,
        state.p_l[lobs_idx].LENGTH_EQ,
        lt=0.0,
        lt_thread=[0.0 for _ in range(state.thread_num)],
    )
    l3 = Line(
        (v1_idx, v2_idx),
        [cobs_idx, cinter2],
        state.p_l[lobs_idx].K1_LENGTH,
        state.p_l[lobs_idx].K1_PCP_LENGTH,
        state.p_l[lobs_idx].K2_LENGTH,
        state.p_l[lobs_idx].LENGTH_EQ,
        lt=0.0,
        lt_thread=[0.0 for _ in range(state.thread_num)],
    )

    a, b = state.p_l[linter2].vi
    state.p_l[linter2].vi = (v1_idx if a == vinter_idx else a, v1_idx if b == vinter_idx else b)
    a, b = state.p_l[linter3].vi
    state.p_l[linter3].vi = (v2_idx if a == vinter_idx else a, v2_idx if b == vinter_idx else b)
    a, b = state.p_l[lobs_idx].vi
    state.p_l[lobs_idx].vi = (v2_idx if a == vobs1 else a, v2_idx if b == vobs1 else b)

    state.p_c[cobs_idx].vi.extend([vinter_idx, v1_idx, v2_idx])
    state.p_c[cobs_idx].li.extend([l1_idx, l2_idx, l3_idx])
    state.p_c[cinter1].vi.append(v1_idx)
    state.p_c[cinter1].li.append(l2_idx)
    _erase_one(state.p_c[cinter2].vi, vinter_idx)
    state.p_c[cinter2].vi.extend([v1_idx, v2_idx])
    state.p_c[cinter2].li.append(l3_idx)

    state.p_v.extend([v1, v2])
    state.p_l.extend([l1, l2, l3])
    if event_log is not None:
        event_log.append(f"{step},intersection,330,{vinter_idx},{cobs_idx},{lobs_idx}")


def cell_intersection(
    state: GlobalState,
    step: int | None = None,
    event_log: list[str] | None = None,
) -> int:
    count = 0
    vlimit = len(state.p_v)
    climit = len(state.p_c)
    for vidx in range(vlimit):
        for cidx in range(climit):
            if len(state.p_v[vidx].ci) == 0:
                break
            if vidx in state.p_c[cidx].vi:
                continue
            if not _contain_vertex_in_cell(state, vidx, cidx):
                continue

            lobs = _closest_line_index(state, vidx, cidx)
            if len(state.p_l[lobs].ci) > 1:
                continue

            vadj = []
            for lidx in state.p_v[vidx].li:
                a, b = state.p_l[lidx].vi
                if a != vidx:
                    vadj.append(a)
                if b != vidx:
                    vadj.append(b)

            deg = len(state.p_v[vidx].li)
            outside = deg
            lobs_touch = 0
            a, b = state.p_l[lobs].vi
            for t in (a, b):
                if t in vadj:
                    outside -= 1
                    lobs_touch += 1
            for vv in state.p_c[cidx].vi:
                if vv in (a, b):
                    continue
                if vv in vadj:
                    outside -= 1
            if deg == 2 and outside == 2 and lobs_touch == 0:
                _intersection_type220(state, vidx, cidx, lobs, step=step, event_log=event_log)
                sort_counter_clockwise(state)
                count += 1
            elif deg == 3 and outside == 2 and lobs_touch == 1:
                _intersection_type321(state, vidx, cidx, lobs, step=step, event_log=event_log)
                sort_counter_clockwise(state)
                count += 1
            elif deg == 2 and outside == 1 and lobs_touch == 1:
                _intersection_type221(state, vidx, cidx, lobs, step=step, event_log=event_log)
                sort_counter_clockwise(state)
                count += 1
            elif deg == 3 and outside == 3 and lobs_touch == 0:
                _intersection_type330(state, vidx, cidx, lobs, step=step, event_log=event_log)
                sort_counter_clockwise(state)
                count += 1
            else:
                import logging

                logging.error(
                    f"Unhandled intersection pattern at step {step}: [{deg}{outside}{lobs_touch}] for vertex {vidx} entering cell {cidx}. deg={deg}, outside={outside}, lobs_touch={lobs_touch}."
                )
    return count
