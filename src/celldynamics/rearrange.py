from __future__ import annotations

from .division import sort_counter_clockwise
from .model import GlobalState


def _erase_one(xs: list[int], x: int) -> None:
    if x in xs:
        xs.remove(x)


def _type442(s: GlobalState, lsep: int) -> None:
    v1, v2 = s.p_l[lsep].vi
    
    ladj1 = ladj2 = -1
    for lidx in s.p_v[v1].li:
        if lidx != lsep and ladj1 == -1:
            ladj1 = lidx
        elif lidx != lsep:
            ladj2 = lidx

    c1 = c2 = c4 = -1
    for i in range(2):
        cidx = s.p_l[ladj1].ci[i] if i < len(s.p_l[ladj1].ci) else -1
        for j in range(2):
            cidx_j = s.p_l[ladj2].ci[j] if j < len(s.p_l[ladj2].ci) else -1
            if cidx != cidx_j:
                continue
            c1 = cidx
            c2 = s.p_l[ladj2].ci[1 - j] if 1 - j < len(s.p_l[ladj2].ci) else -1
            c4 = s.p_l[ladj1].ci[1 - i] if 1 - i < len(s.p_l[ladj1].ci) else -1

    c3 = -1
    for cidx in s.p_v[v2].ci:
        if cidx != c2 and cidx != c4:
            c3 = cidx

    ladj3 = -1
    for lidx in s.p_v[v2].li:
        if lidx == lsep:
            continue
        ci = s.p_l[lidx].ci
        if (len(ci) > 0 and ci[0] == c2) or (len(ci) > 1 and ci[1] == c2):
            ladj3 = lidx

    p1 = s.p_v[v1].loc[0]
    p2 = s.p_v[v2].loc[0]
    center = (p1 + p2) * 0.5
    d1 = p1 - center
    d2 = p2 - center
    
    vec_v = p2 - p1
    # Check orientation using ladj1's other vertex
    v_adj1 = s.p_l[ladj1].vi[0] if s.p_l[ladj1].vi[1] == v1 else s.p_l[ladj1].vi[1]
    vec_adj1 = s.p_v[v_adj1].loc[0] - p1
    cross_val = vec_v.x * vec_adj1.y - vec_v.y * vec_adj1.x
    rot_sign = 1.0 if cross_val > 0 else -1.0

    _erase_one(s.p_v[v1].li, ladj1)
    s.p_v[v1].li.append(ladj3)
    _erase_one(s.p_v[v1].ci, c4)
    s.p_v[v1].ci.append(c3)
    s.p_v[v1].loc[0].x = center.x - rot_sign * d1.y
    s.p_v[v1].loc[0].y = center.y + rot_sign * d1.x

    _erase_one(s.p_v[v2].li, ladj3)
    s.p_v[v2].li.append(ladj1)
    _erase_one(s.p_v[v2].ci, c2)
    s.p_v[v2].ci.append(c1)
    s.p_v[v2].loc[0].x = center.x - rot_sign * d2.y
    s.p_v[v2].loc[0].y = center.y + rot_sign * d2.x

    _erase_one(s.p_l[lsep].ci, c2)
    _erase_one(s.p_l[lsep].ci, c4)
    s.p_l[lsep].ci.append(c1)
    s.p_l[lsep].ci.append(c3)

    a, b = s.p_l[ladj1].vi
    s.p_l[ladj1].vi = (v2 if a == v1 else a, v2 if b == v1 else b)
    a, b = s.p_l[ladj3].vi
    s.p_l[ladj3].vi = (v1 if a == v2 else a, v1 if b == v2 else b)

    s.p_c[c1].vi.append(v2)
    s.p_c[c1].li.append(lsep)
    _erase_one(s.p_c[c2].vi, v2)
    _erase_one(s.p_c[c2].li, lsep)
    s.p_c[c3].vi.append(v1)
    s.p_c[c3].li.append(lsep)
    _erase_one(s.p_c[c4].vi, v1)
    _erase_one(s.p_c[c4].li, lsep)


def _type432(s: GlobalState, lsep: int) -> None:
    v1, v2 = s.p_l[lsep].vi
    
    ladj1 = ladj2 = -1
    for lidx in s.p_v[v1].li:
        if lidx != lsep and ladj1 == -1:
            ladj1 = lidx
        elif lidx != lsep:
            ladj2 = lidx

    c1 = c2 = c4 = -1
    for i in range(2):
        cidx = s.p_l[ladj1].ci[i] if i < len(s.p_l[ladj1].ci) else -1
        for j in range(2):
            cidx_j = s.p_l[ladj2].ci[j] if j < len(s.p_l[ladj2].ci) else -1
            if cidx != cidx_j:
                continue
            c1 = cidx
            c2 = s.p_l[ladj2].ci[1 - j] if 1 - j < len(s.p_l[ladj2].ci) else -1
            c4 = s.p_l[ladj1].ci[1 - i] if 1 - i < len(s.p_l[ladj1].ci) else -1

    ladj3 = -1
    for lidx in s.p_v[v2].li:
        if lidx != lsep and len(s.p_l[lidx].ci) > 0 and s.p_l[lidx].ci[0] == c2:
            ladj3 = lidx

    p1 = s.p_v[v1].loc[0]
    p2 = s.p_v[v2].loc[0]
    center = (p1 + p2) * 0.5
    d1 = p1 - center
    d2 = p2 - center
    
    vec_v = p2 - p1
    # Check orientation using ladj1's other vertex
    v_adj1 = s.p_l[ladj1].vi[0] if s.p_l[ladj1].vi[1] == v1 else s.p_l[ladj1].vi[1]
    vec_adj1 = s.p_v[v_adj1].loc[0] - p1
    cross_val = vec_v.x * vec_adj1.y - vec_v.y * vec_adj1.x
    rot_sign = 1.0 if cross_val > 0 else -1.0

    _erase_one(s.p_v[v1].li, ladj1)
    s.p_v[v1].li.append(ladj3)
    _erase_one(s.p_v[v1].ci, c4)
    s.p_v[v1].loc[0].x = center.x - rot_sign * d1.y
    s.p_v[v1].loc[0].y = center.y + rot_sign * d1.x

    _erase_one(s.p_v[v2].li, ladj3)
    s.p_v[v2].li.append(ladj1)
    _erase_one(s.p_v[v2].ci, c2)
    s.p_v[v2].ci.append(c1)
    s.p_v[v2].loc[0].x = center.x - rot_sign * d2.y
    s.p_v[v2].loc[0].y = center.y + rot_sign * d2.x

    _erase_one(s.p_l[lsep].ci, c2)
    _erase_one(s.p_l[lsep].ci, c4)
    s.p_l[lsep].ci.append(c1)

    a, b = s.p_l[ladj1].vi
    s.p_l[ladj1].vi = (v2 if a == v1 else a, v2 if b == v1 else b)
    a, b = s.p_l[ladj3].vi
    s.p_l[ladj3].vi = (v1 if a == v2 else a, v1 if b == v2 else b)

    s.p_c[c1].vi.append(v2)
    s.p_c[c1].li.append(lsep)
    _erase_one(s.p_c[c2].vi, v2)
    _erase_one(s.p_c[c2].li, lsep)
    _erase_one(s.p_c[c4].vi, v1)
    _erase_one(s.p_c[c4].li, lsep)


def _type431(s: GlobalState, lsep: int) -> None:
    v1, v2 = s.p_l[lsep].vi
    
    ladj1 = ladj2 = -1
    for lidx in s.p_v[v1].li:
        if lidx != lsep and len(s.p_l[lidx].ci) == 2:
            ladj1 = lidx
        elif lidx != lsep:
            ladj2 = lidx

    c1 = s.p_l[ladj2].ci[0]
    c4 = s.p_l[lsep].ci[0]

    c3 = -1
    for cidx in s.p_v[v2].ci:
        if cidx != c4:
            c3 = cidx

    ladj3 = -1
    for lidx in s.p_v[v2].li:
        if lidx != lsep and len(s.p_l[lidx].ci) == 1:
            ladj3 = lidx

    p1 = s.p_v[v1].loc[0]
    p2 = s.p_v[v2].loc[0]
    center = (p1 + p2) * 0.5
    d1 = p1 - center
    d2 = p2 - center
    
    vec_v = p2 - p1
    # Check orientation using ladj1's other vertex
    v_adj1 = s.p_l[ladj1].vi[0] if s.p_l[ladj1].vi[1] == v1 else s.p_l[ladj1].vi[1]
    vec_adj1 = s.p_v[v_adj1].loc[0] - p1
    cross_val = vec_v.x * vec_adj1.y - vec_v.y * vec_adj1.x
    rot_sign = 1.0 if cross_val > 0 else -1.0

    _erase_one(s.p_v[v1].li, ladj1)
    s.p_v[v1].li.append(ladj3)
    _erase_one(s.p_v[v1].ci, c4)
    s.p_v[v1].ci.append(c3)
    s.p_v[v1].loc[0].x = center.x - rot_sign * d1.y
    s.p_v[v1].loc[0].y = center.y + rot_sign * d1.x

    _erase_one(s.p_v[v2].li, ladj3)
    s.p_v[v2].li.append(ladj1)
    s.p_v[v2].ci.append(c1)
    s.p_v[v2].loc[0].x = center.x - rot_sign * d2.y
    s.p_v[v2].loc[0].y = center.y + rot_sign * d2.x

    _erase_one(s.p_l[lsep].ci, c4)
    s.p_l[lsep].ci.append(c1)
    s.p_l[lsep].ci.append(c3)

    a, b = s.p_l[ladj1].vi
    s.p_l[ladj1].vi = (v2 if a == v1 else a, v2 if b == v1 else b)
    a, b = s.p_l[ladj3].vi
    s.p_l[ladj3].vi = (v1 if a == v2 else a, v1 if b == v2 else b)

    s.p_c[c1].vi.append(v2)
    s.p_c[c1].li.append(lsep)
    s.p_c[c3].vi.append(v1)
    s.p_c[c3].li.append(lsep)
    _erase_one(s.p_c[c4].vi, v1)
    _erase_one(s.p_c[c4].li, lsep)


def _type332(s: GlobalState, lsep: int) -> None:
    v1, v2 = s.p_l[lsep].vi
    
    ladj1 = ladj2 = -1
    for lidx in s.p_v[v1].li:
        if lidx != lsep and ladj1 == -1:
            ladj1 = lidx
        elif lidx != lsep:
            ladj2 = lidx

    c1 = c2 = c3 = -1
    for i in range(2):
        cidx = s.p_l[ladj1].ci[i] if i < len(s.p_l[ladj1].ci) else -1
        for j in range(2):
            cidx_j = s.p_l[ladj2].ci[j] if j < len(s.p_l[ladj2].ci) else -1
            if cidx != cidx_j:
                continue
            c1 = cidx
            c2 = s.p_l[ladj2].ci[1 - j] if 1 - j < len(s.p_l[ladj2].ci) else -1
            c3 = s.p_l[ladj1].ci[1 - i] if 1 - i < len(s.p_l[ladj1].ci) else -1

    ladj3 = -1
    for lidx in s.p_v[v2].li:
        if lidx != lsep:
            ladj3 = lidx

    p1 = s.p_v[v1].loc[0]
    p2 = s.p_v[v2].loc[0]
    center = (p1 + p2) * 0.5
    d1 = p1 - center
    d2 = p2 - center
    
    vec_v = p2 - p1
    # Check orientation using ladj1's other vertex
    v_adj1 = s.p_l[ladj1].vi[0] if s.p_l[ladj1].vi[1] == v1 else s.p_l[ladj1].vi[1]
    vec_adj1 = s.p_v[v_adj1].loc[0] - p1
    cross_val = vec_v.x * vec_adj1.y - vec_v.y * vec_adj1.x
    rot_sign = 1.0 if cross_val > 0 else -1.0

    _erase_one(s.p_v[v1].li, ladj1)
    s.p_v[v1].li.append(ladj3)
    s.p_v[v1].loc[0].x = center.x - rot_sign * d1.y
    s.p_v[v1].loc[0].y = center.y + rot_sign * d1.x

    _erase_one(s.p_v[v2].li, ladj3)
    s.p_v[v2].li.append(ladj1)
    _erase_one(s.p_v[v2].ci, c2)
    s.p_v[v2].ci.append(c1)
    s.p_v[v2].loc[0].x = center.x - rot_sign * d2.y
    s.p_v[v2].loc[0].y = center.y + rot_sign * d2.x

    _erase_one(s.p_l[lsep].ci, c2)
    s.p_l[lsep].ci.append(c1)

    a, b = s.p_l[ladj1].vi
    s.p_l[ladj1].vi = (v2 if a == v1 else a, v2 if b == v1 else b)
    a, b = s.p_l[ladj3].vi
    s.p_l[ladj3].vi = (v1 if a == v2 else a, v1 if b == v2 else b)

    s.p_c[c1].vi.append(v2)
    s.p_c[c1].li.append(lsep)
    _erase_one(s.p_c[c2].vi, v2)
    _erase_one(s.p_c[c2].li, lsep)


def _type321(s: GlobalState, lsep: int) -> None:
    v1, v2 = s.p_l[lsep].vi
    
    ladj1 = ladj2 = -1
    for lidx in s.p_v[v1].li:
        if len(s.p_l[lidx].ci) == 2:
            ladj1 = lidx
        elif lidx != lsep:
            ladj2 = lidx

    c1 = s.p_l[ladj2].ci[0]
    c4 = s.p_l[lsep].ci[0]

    p1 = s.p_v[v1].loc[0]
    p2 = s.p_v[v2].loc[0]
    center = (p1 + p2) * 0.5
    d1 = p1 - center
    d2 = p2 - center
    
    vec_v = p2 - p1
    # Check orientation using ladj1's other vertex
    v_adj1 = s.p_l[ladj1].vi[0] if s.p_l[ladj1].vi[1] == v1 else s.p_l[ladj1].vi[1]
    vec_adj1 = s.p_v[v_adj1].loc[0] - p1
    cross_val = vec_v.x * vec_adj1.y - vec_v.y * vec_adj1.x
    rot_sign = 1.0 if cross_val > 0 else -1.0

    _erase_one(s.p_v[v1].li, ladj1)
    _erase_one(s.p_v[v1].ci, c4)
    s.p_v[v1].loc[0].x = center.x - rot_sign * d1.y
    s.p_v[v1].loc[0].y = center.y + rot_sign * d1.x

    s.p_v[v2].li.append(ladj1)
    s.p_v[v2].ci.append(c1)
    s.p_v[v2].loc[0].x = center.x - rot_sign * d2.y
    s.p_v[v2].loc[0].y = center.y + rot_sign * d2.x

    _erase_one(s.p_l[lsep].ci, c4)
    s.p_l[lsep].ci.append(c1)

    a, b = s.p_l[ladj1].vi
    s.p_l[ladj1].vi = (v2 if a == v1 else a, v2 if b == v1 else b)

    s.p_c[c1].vi.append(v2)
    s.p_c[c1].li.append(lsep)
    _erase_one(s.p_c[c4].vi, v1)
    _erase_one(s.p_c[c4].li, lsep)


def _type211(s: GlobalState, lsep: int) -> None:
    v1, v2 = s.p_l[lsep].vi
    
    ladj2 = -1
    for lidx in s.p_v[v2].li:
        if lidx != lsep:
            ladj2 = lidx

    cadj = s.p_l[lsep].ci[0]

    p1 = s.p_v[v1].loc[0]
    p2 = s.p_v[v2].loc[0]
    center = (p1 + p2) * 0.5

    _erase_one(s.p_v[v1].li, lsep)
    s.p_v[v1].li.append(ladj2)
    s.p_v[v1].loc[0] = center

    s.p_v[v2].li.clear()
    s.p_v[v2].ci.clear()

    s.p_l[lsep].ci.clear()

    a, b = s.p_l[ladj2].vi
    s.p_l[ladj2].vi = (v1 if a == v2 else a, v1 if b == v2 else b)

    _erase_one(s.p_c[cadj].vi, v2)
    _erase_one(s.p_c[cadj].li, lsep)



def cell_rearrange2(
    state: GlobalState,
    threshold: float,
    step: int | None = None,
    event_log: list[str] | None = None,
    vacant_area_th: float | None = None,
) -> int:
    count = 0
    if vacant_area_th is not None:
        c_vac = _detect_and_fix_triangle_vacants(state, vacant_area_th, step, event_log)
        if c_vac > 0:
            sort_counter_clockwise(state)
        count += c_vac


    for lsep in range(len(state.p_l)):
        lp = state.p_l[lsep]
        v1, v2 = lp.vi
        if len(lp.ci) == 0:
            continue
        if (state.p_v[v1].loc[0] - state.p_v[v2].loc[0]).norm() > threshold:
            continue
        if any(len(state.p_c[c].vi) <= 3 for c in lp.ci):
            continue

        if len(state.p_v[v1].li) < len(state.p_v[v2].li):
            lp.vi = (v2, v1)
            v1, v2 = v2, v1
        elif len(state.p_v[v1].li) == len(state.p_v[v2].li) and len(state.p_v[v1].ci) < len(state.p_v[v2].ci):
            lp.vi = (v2, v1)
            v1, v2 = v2, v1

        line_cnt = len(state.p_v[v1].li) + len(state.p_v[v2].li) - 2
        cell_total = len(state.p_v[v1].ci) + len(state.p_v[v2].ci) - len(lp.ci)
        cell_lsep = len(lp.ci)

        changed = False
        if line_cnt == 4 and cell_total == 4 and cell_lsep == 2:
            _type442(state, lsep)
            count += 1
            changed = True
            if event_log is not None:
                event_log.append(f"{step},rearrange,442,{lsep}")
        elif line_cnt == 4 and cell_total == 3 and cell_lsep == 2:
            _type432(state, lsep)
            count += 1
            changed = True
            if event_log is not None:
                event_log.append(f"{step},rearrange,432,{lsep}")
        elif line_cnt == 4 and cell_total == 3 and cell_lsep == 1:
            _type431(state, lsep)
            count += 1
            changed = True
            if event_log is not None:
                event_log.append(f"{step},rearrange,431,{lsep}")
        elif line_cnt == 3 and cell_total == 3 and cell_lsep == 2:
            _type332(state, lsep)
            count += 1
            changed = True
            if event_log is not None:
                event_log.append(f"{step},rearrange,332,{lsep}")
        elif line_cnt == 3 and cell_total == 2 and cell_lsep == 1:
            _type321(state, lsep)
            count += 1
            changed = True
            if event_log is not None:
                event_log.append(f"{step},rearrange,321,{lsep}")
        elif line_cnt == 2 and cell_total == 1 and cell_lsep == 1:
            _type211(state, lsep)
            count += 1
            changed = True
            if event_log is not None:
                event_log.append(f"{step},rearrange,211,{lsep}")

        if changed:
            sort_counter_clockwise(state)
    return count
def _detect_and_fix_triangle_vacants(
    state: GlobalState, a_th: float, step: int | None, event_log: list[str] | None
) -> int:
    count = 0
    active_bnd_lines = set()
    for lidx, line in enumerate(state.p_l):
        if len(line.ci) == 1 and len(line.vi) == 2 and line.vi[0] != line.vi[1]:
            active_bnd_lines.add(lidx)

    bnd_adj: dict[int, list[int]] = {}
    for lidx in active_bnd_lines:
        v1, v2 = state.p_l[lidx].vi
        bnd_adj.setdefault(v1, []).append(lidx)
        bnd_adj.setdefault(v2, []).append(lidx)

    triangles: list[tuple[int, int, int]] = []
    visited_triangles = set()

    for v1, lines1 in bnd_adj.items():
        if len(lines1) < 2:
            continue
        for i in range(len(lines1)):
            for j in range(i + 1, len(lines1)):
                l1 = lines1[i]
                l2 = lines1[j]
                
                v2 = state.p_l[l1].vi[0] if state.p_l[l1].vi[1] == v1 else state.p_l[l1].vi[1]
                v3 = state.p_l[l2].vi[0] if state.p_l[l2].vi[1] == v1 else state.p_l[l2].vi[1]
                
                if v2 == v3:
                    continue
                
                common_lines = set(bnd_adj.get(v2, [])) & set(bnd_adj.get(v3, []))
                for l3 in common_lines:
                    if l3 in active_bnd_lines:
                        tri = tuple(sorted([l1, l2, l3]))
                        if tri not in visited_triangles:
                            visited_triangles.add(tri)
                            triangles.append(tri)

    deactivated_lines = set()
    deactivated_vertices = set()

    for tri_lines in triangles:
        l1, l2, l3 = tri_lines
        if l1 in deactivated_lines or l2 in deactivated_lines or l3 in deactivated_lines:
            continue
            
        v_set = set(state.p_l[l1].vi) | set(state.p_l[l2].vi) | set(state.p_l[l3].vi)
        if len(v_set) != 3:
            continue
        
        v_list = list(v_set)
        vA, vB, vC = v_list[0], v_list[1], v_list[2]
        
        pA = state.p_v[vA].loc[0]
        pB = state.p_v[vB].loc[0]
        pC = state.p_v[vC].loc[0]
        
        vec1 = pB - pA
        vec2 = pC - pA
        area = 0.5 * abs(vec1.x * vec2.y - vec1.y * vec2.x)
        
        if area > a_th:
            continue
            
        count += 1
        
        G = (pA + pB + pC) * (1.0 / 3.0)
        
        v_keep = vA
        v_drop = {vB, vC}
        
        state.p_v[v_keep].loc[0] = G
        
        lines_to_rewire = set()
        for v in v_set:
            for lidx in state.p_v[v].li:
                if lidx not in tri_lines:
                    lines_to_rewire.add(lidx)
                    
        for lidx in lines_to_rewire:
            vi = list(state.p_l[lidx].vi)
            changed_l = False
            for i in range(2):
                if vi[i] in v_drop:
                    vi[i] = v_keep
                    changed_l = True
            if changed_l:
                state.p_l[lidx].vi = (vi[0], vi[1])
        
        state.p_v[v_keep].li = list(lines_to_rewire)
        
        c_set = set()
        for v in v_set:
            c_set.update(state.p_v[v].ci)
            
        state.p_v[v_keep].ci = list(c_set)
        
        for cidx in c_set:
            cell = state.p_c[cidx]
            
            new_vi = []
            for v in cell.vi:
                nv = v_keep if v in v_drop else v
                if not new_vi or new_vi[-1] != nv:
                    new_vi.append(nv)
            if len(new_vi) > 1 and new_vi[0] == new_vi[-1]:
                new_vi.pop()
            
            seen = set()
            final_vi = []
            for nv in new_vi:
                if nv not in seen:
                    seen.add(nv)
                    final_vi.append(nv)
            
            cell.vi = final_vi
            cell.li = [lidx for lidx in cell.li if lidx not in tri_lines]
            
        for lidx in tri_lines:
            state.p_l[lidx].ci.clear()
            deactivated_lines.add(lidx)
            
        for v in v_drop:
            state.p_v[v].li.clear()
            state.p_v[v].ci.clear()
            deactivated_vertices.add(v)
            
        if event_log is not None:
            event_log.append(f"{step},vacant_triangle,{area:.6f},{l1},{l2},{l3}")

    return count
