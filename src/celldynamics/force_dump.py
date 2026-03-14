from __future__ import annotations

from pathlib import Path

from .model import GlobalState


def _zpad(step: int, width: int = 10) -> str:
    return f"{step:0{width}d}"


def dump_force_snapshot(
    state: GlobalState,
    deg: int,
    phase: str,
    step: int,
    out_dir: str | Path,
) -> tuple[Path, Path]:
    """Write C++-compatible force/lt CSV snapshots for parity comparison."""
    d = Path(out_dir)
    d.mkdir(parents=True, exist_ok=True)

    s = _zpad(step)
    f_force = d / f"force_{phase}_{s}.csv"
    f_lt = d / f"lt_{phase}_{s}.csv"

    with f_force.open("w", encoding="utf-8") as ofs:
        ofs.write("step,deg,vertex_id,fx,fy,fz\n")
        for i, vp in enumerate(state.p_v):
            ofs.write(
                f"{step},{deg},{i},"
                f"{vp.frc[deg].x:.17g},{vp.frc[deg].y:.17g},{vp.frc[deg].z:.17g}\n"
            )

    with f_lt.open("w", encoding="utf-8") as ofs:
        ofs.write("step,line_id,lt\n")
        for i, lp in enumerate(state.p_l):
            ofs.write(f"{step},{i},{lp.lt:.17g}\n")

    return f_force, f_lt
