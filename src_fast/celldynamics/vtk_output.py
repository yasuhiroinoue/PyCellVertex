from __future__ import annotations

import math
from pathlib import Path

from .model import GlobalState


def _cxx_float(v: float) -> str:
    """Approximate default C++ ostream float formatting (precision=6, defaultfloat)."""
    return format(v, ".6g")


def _header() -> list[str]:
    return [
        "# vtk DataFile Version 2.0",
        "2D-vertex",
        "ASCII",
        "DATASET UNSTRUCTURED_GRID",
    ]


def render_line_vtk(state: GlobalState, alive_only: bool = False) -> str:
    lines = _header()
    
    if alive_only:
        alive_lines = [lp for lp in state.p_l if len(lp.ci) > 0]
        alive_vids = set()
        for lp in alive_lines:
            alive_vids.add(lp.vi[0])
            alive_vids.add(lp.vi[1])
        
        sorted_vids = sorted(list(alive_vids))
        vid_map = {old: new for new, old in enumerate(sorted_vids)}
        
        lines.append(f"POINTS {len(sorted_vids)} float")
        for vid in sorted_vids:
            vp = state.p_v[vid]
            lines.append(f"{_cxx_float(vp.loc[0].x)} {_cxx_float(vp.loc[0].y)} {_cxx_float(vp.loc[0].z)}")
            
        lines.append(f"CELLS {len(alive_lines)} {len(alive_lines) * 3}")
        for lp in alive_lines:
            lines.append(f"2 {vid_map[lp.vi[0]]} {vid_map[lp.vi[1]]}")
            
        lines.append(f"CELL_TYPES {len(alive_lines)}")
        for _ in alive_lines:
            lines.append("3")
            
        lines.append(f"CELL_DATA {len(alive_lines)}")
        lines.append("SCALARS line_tension float")
        lines.append("LOOKUP_TABLE default")
        for lp in alive_lines:
            lines.append(_cxx_float(lp.lt))
    else:
        lines.append(f"POINTS {len(state.p_v)} float")
        for vp in state.p_v:
            lines.append(f"{_cxx_float(vp.loc[0].x)} {_cxx_float(vp.loc[0].y)} {_cxx_float(vp.loc[0].z)}")

        lines.append(f"CELLS {len(state.p_l)} {len(state.p_l) * 3}")
        for lp in state.p_l:
            lines.append(f"2 {lp.vi[0]} {lp.vi[1]}")

        lines.append(f"CELL_TYPES {len(state.p_l)}")
        for _ in state.p_l:
            lines.append("3")

        lines.append(f"CELL_DATA {len(state.p_l)}")
        lines.append("SCALARS line_tension float")
        lines.append("LOOKUP_TABLE default")
        for lp in state.p_l:
            lines.append(_cxx_float(lp.lt))

    return "\n".join(lines) + "\n"


def render_face_vtk(state: GlobalState, alive_only: bool = False) -> str:
    lines = _header()
    
    if alive_only:
        alive_cells = state.p_c  # as per instructions, use all cells in state
        alive_vids = set()
        for cp in alive_cells:
            for vid in cp.vi:
                alive_vids.add(vid)
                
        sorted_vids = sorted(list(alive_vids))
        vid_map = {old: new for new, old in enumerate(sorted_vids)}
        
        lines.append(f"POINTS {len(sorted_vids)} float")
        for vid in sorted_vids:
            vp = state.p_v[vid]
            lines.append(f"{_cxx_float(vp.loc[0].x)} {_cxx_float(vp.loc[0].y)} {_cxx_float(vp.loc[0].z)}")
            
        cells_size = sum(len(cp.vi) + 1 for cp in alive_cells)
        lines.append(f"CELLS {len(alive_cells)} {cells_size}")
        for cp in alive_cells:
            if cp.vi:
                lines.append(f"{len(cp.vi)} {' '.join(str(vid_map[v]) for v in cp.vi)}")
            else:
                lines.append("0")

        lines.append(f"CELL_TYPES {len(alive_cells)}")
        for _ in alive_cells:
            lines.append("7")

        lines.append(f"CELL_DATA {len(alive_cells)}")
        lines.append("SCALARS phase float")
        lines.append("LOOKUP_TABLE default")
        for cp in alive_cells:
            sin_t = math.sin(2.0 * math.pi * cp.cell_time / cp.cell_T - cp.cell_phase)
            lines.append(_cxx_float(sin_t))
    else:
        lines.append(f"POINTS {len(state.p_v)} float")
        for vp in state.p_v:
            lines.append(f"{_cxx_float(vp.loc[0].x)} {_cxx_float(vp.loc[0].y)} {_cxx_float(vp.loc[0].z)}")

        cells_size = sum(len(cp.vi) + 1 for cp in state.p_c)
        lines.append(f"CELLS {len(state.p_c)} {cells_size}")
        for cp in state.p_c:
            if cp.vi:
                lines.append(f"{len(cp.vi)} {' '.join(str(v) for v in cp.vi)}")
            else:
                lines.append("0")

        lines.append(f"CELL_TYPES {len(state.p_c)}")
        for _ in state.p_c:
            lines.append("7")

        lines.append(f"CELL_DATA {len(state.p_c)}")
        lines.append("SCALARS phase float")
        lines.append("LOOKUP_TABLE default")
        for cp in state.p_c:
            sin_t = math.sin(2.0 * math.pi * cp.cell_time / cp.cell_T - cp.cell_phase)
            lines.append(_cxx_float(sin_t))

    return "\n".join(lines) + "\n"


def dump_vtk_snapshot(state: GlobalState, step: int, out_dir: Path, step_scale: int = 10000, alive_only: bool = False, vtk_mode: str = "raw") -> tuple[Path, Path]:
    if alive_only and vtk_mode == "raw":
        vtk_mode = "alive"

    token = f"{step:010d}"
    last_line = None
    last_face = None

    if vtk_mode in ("raw", "both"):
        raw_dir = out_dir / "vtk_raw"
        raw_dir.mkdir(parents=True, exist_ok=True)
        f_line = raw_dir / f"2dv_line{token}.vtk"
        f_face = raw_dir / f"2dv_face{token}.vtk"
        f_line.write_text(render_line_vtk(state, alive_only=False), encoding="utf-8")
        f_face.write_text(render_face_vtk(state, alive_only=False), encoding="utf-8")
        last_line, last_face = f_line, f_face

    if vtk_mode in ("alive", "both"):
        alive_dir = out_dir / "vtk_alive"
        alive_dir.mkdir(parents=True, exist_ok=True)
        f_line = alive_dir / f"2dv_line{token}.vtk"
        f_face = alive_dir / f"2dv_face{token}.vtk"
        f_line.write_text(render_line_vtk(state, alive_only=True), encoding="utf-8")
        f_face.write_text(render_face_vtk(state, alive_only=True), encoding="utf-8")
        last_line, last_face = f_line, f_face

    return last_line, last_face
