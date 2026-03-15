from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv

import numpy as np

from .validation import Tolerance


@dataclass(frozen=True)
class ForceMismatch:
    filename: str
    row: int
    col: str
    actual: float
    expected: float
    abs_diff: float


def _read_force_csv(path: Path) -> tuple[list[str], np.ndarray]:
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        cols = ["fx", "fy", "fz"] if "fx" in reader.fieldnames else ["lt"]
        rows = []
        for row in reader:
            rows.append([float(row[c]) for c in cols])
    return cols, np.asarray(rows, dtype=float)


def first_force_mismatch(
    expected_dir: Path,
    actual_dir: Path,
    tol: Tolerance = Tolerance(),
    prefix: str = "force_",
    allow_missing: bool = False,
) -> ForceMismatch | None:
    exp_files = sorted([p for p in expected_dir.iterdir() if p.is_file() and p.name.startswith(prefix)])
    act_map = {p.name: p for p in actual_dir.iterdir() if p.is_file() and p.name.startswith(prefix)}

    for ef in exp_files:
        af = act_map.get(ef.name)
        if af is None:
            if allow_missing:
                continue
            return ForceMismatch(ef.name, -1, "missing", np.nan, np.nan, np.inf)

        cols, e = _read_force_csv(ef)
        cols2, a = _read_force_csv(af)
        if cols != cols2 or e.shape != a.shape:
            return ForceMismatch(ef.name, -2, "shape", np.nan, np.nan, np.inf)

        if not np.allclose(a, e, rtol=tol.rtol, atol=tol.atol):
            d = np.abs(a - e)
            i = np.unravel_index(np.argmax(d), d.shape)
            return ForceMismatch(
                filename=ef.name,
                row=int(i[0]),
                col=cols[int(i[1])],
                actual=float(a[i]),
                expected=float(e[i]),
                abs_diff=float(d[i]),
            )

    return None
