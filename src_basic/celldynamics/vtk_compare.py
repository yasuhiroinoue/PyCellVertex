from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .validation import Tolerance


_STEP_RE = re.compile(r"(\d+)")


@dataclass(frozen=True)
class NumericMismatch:
    filename: str
    max_abs_diff: float
    index: tuple[int, int]


@dataclass(frozen=True)
class VtkStringMismatch:
    filename: str
    kind: str
    step: int
    reason: str
    left: str | None
    right: str | None


@dataclass(frozen=True)
class VtkNumericMismatch:
    filename: str
    kind: str
    step: int
    section: str
    reason: str
    index: tuple[int, ...] | None
    left_value: float | int | str | None
    right_value: float | int | str | None
    max_abs_diff: float | None = None


@dataclass(frozen=True)
class _ParsedVtk:
    points: np.ndarray
    cells: np.ndarray
    cell_types: np.ndarray
    cell_data: dict[str, np.ndarray]


def _extract_step(name: str) -> int:
    m = _STEP_RE.findall(name)
    return int(m[-1]) if m else -1


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _collect_numeric(lines: list[str], start: int, expected: int, cast: type[float] | type[int]) -> tuple[np.ndarray, int]:
    vals: list[float | int] = []
    i = start
    while i < len(lines) and len(vals) < expected:
        s = lines[i].strip()
        if not s:
            i += 1
            continue
        vals.extend(cast(x) for x in s.split())
        i += 1
    if len(vals) < expected:
        raise ValueError(f"expected {expected} values, got {len(vals)}")
    return np.asarray(vals[:expected]), i


def _parse_vtk(path: Path) -> _ParsedVtk:
    lines = path.read_text(encoding="utf-8").splitlines()

    points = None
    cells = None
    cell_types = None
    cell_data: dict[str, np.ndarray] = {}
    n_cell_data: int | None = None

    i = 0
    while i < len(lines):
        s = lines[i].strip()
        if not s:
            i += 1
            continue

        if s.startswith("POINTS "):
            parts = s.split()
            n_points = int(parts[1])
            arr, i = _collect_numeric(lines, i + 1, n_points * 3, float)
            points = arr.reshape(n_points, 3)
            continue

        if s.startswith("CELLS "):
            parts = s.split()
            total_size = int(parts[2])
            arr, i = _collect_numeric(lines, i + 1, total_size, int)
            cells = arr.astype(int)
            continue

        if s.startswith("CELL_TYPES "):
            parts = s.split()
            n_types = int(parts[1])
            arr, i = _collect_numeric(lines, i + 1, n_types, int)
            cell_types = arr.astype(int)
            continue

        if s.startswith("CELL_DATA "):
            n_cell_data = int(s.split()[1])
            i += 1
            continue

        if s.startswith("SCALARS "):
            if n_cell_data is None:
                raise ValueError("SCALARS appeared before CELL_DATA")
            parts = s.split()
            name = parts[1]
            i += 1
            if i >= len(lines) or not lines[i].strip().startswith("LOOKUP_TABLE"):
                raise ValueError(f"SCALARS {name} missing LOOKUP_TABLE")
            arr, i = _collect_numeric(lines, i + 1, n_cell_data, float)
            cell_data[name] = arr.astype(float)
            continue

        i += 1

    if points is None:
        raise ValueError(f"POINTS section not found: {path}")
    if cells is None:
        raise ValueError(f"CELLS section not found: {path}")
    if cell_types is None:
        raise ValueError(f"CELL_TYPES section not found: {path}")

    return _ParsedVtk(points=points, cells=cells, cell_types=cell_types, cell_data=cell_data)


def read_vtk_points(path: Path) -> np.ndarray:
    return _parse_vtk(path).points


def _first_allclose_mismatch(a: np.ndarray, b: np.ndarray, tol: Tolerance) -> tuple[tuple[int, ...], float, float, float] | None:
    if a.shape != b.shape:
        return None
    if np.allclose(a, b, rtol=tol.rtol, atol=tol.atol):
        return None
    diff = np.abs(a - b)
    idx = np.unravel_index(np.argmax(diff), diff.shape)
    return tuple(int(i) for i in idx), float(diff[idx]), float(a[idx]), float(b[idx])


def _first_exact_mismatch(a: np.ndarray, b: np.ndarray) -> tuple[int, int, int] | None:
    if a.shape != b.shape:
        return None
    neq = np.flatnonzero(a != b)
    if len(neq) == 0:
        return None
    i = int(neq[0])
    return i, int(a[i]), int(b[i])


def first_numeric_mismatch_by_step(
    left_dir: Path,
    right_dir: Path,
    tol: Tolerance = Tolerance(),
    prefix: str = "2dv_face",
) -> NumericMismatch | None:
    left = sorted([p for p in left_dir.iterdir() if p.is_file() and p.name.startswith(prefix)])
    right_map = {p.name: p for p in right_dir.iterdir() if p.is_file() and p.name.startswith(prefix)}

    for lp in left:
        rp = right_map.get(lp.name)
        if rp is None:
            return NumericMismatch(lp.name, float("inf"), (-1, -1))

        a = read_vtk_points(lp)
        b = read_vtk_points(rp)
        if a.shape != b.shape:
            return NumericMismatch(lp.name, float("inf"), (-2, -2))

        if not np.allclose(a, b, rtol=tol.rtol, atol=tol.atol):
            diff = np.abs(a - b)
            i = np.unravel_index(np.argmax(diff), diff.shape)
            return NumericMismatch(lp.name, float(diff[i]), (int(i[0]), int(i[1])))

    return None


def first_vtk_string_mismatch_by_step(left_dir: Path, right_dir: Path, kind: str = "face") -> VtkStringMismatch | None:
    if kind not in {"face", "line"}:
        raise ValueError(f"kind must be 'face' or 'line': {kind}")

    prefix = f"2dv_{kind}"
    left_files = {p.name: p for p in left_dir.iterdir() if p.is_file() and p.name.startswith(prefix) and p.suffix == ".vtk"}
    right_files = {p.name: p for p in right_dir.iterdir() if p.is_file() and p.name.startswith(prefix) and p.suffix == ".vtk"}

    all_names = sorted(set(left_files) | set(right_files), key=lambda n: (_extract_step(n), n))
    for name in all_names:
        l = left_files.get(name)
        r = right_files.get(name)
        step = _extract_step(name)
        if l is None:
            return VtkStringMismatch(name, kind, step, "missing_in_left", None, str(r))
        if r is None:
            return VtkStringMismatch(name, kind, step, "missing_in_right", str(l), None)

        hl = _sha256_file(l)
        hr = _sha256_file(r)
        if hl != hr:
            return VtkStringMismatch(name, kind, step, "content_diff", hl, hr)

    return None


def first_vtk_numeric_mismatch_by_step(
    left_dir: Path,
    right_dir: Path,
    kind: str = "face",
    tol: Tolerance = Tolerance(),
    check_topology: bool = True,
) -> VtkNumericMismatch | None:
    if kind not in {"face", "line"}:
        raise ValueError(f"kind must be 'face' or 'line': {kind}")

    prefix = f"2dv_{kind}"
    left_files = {p.name: p for p in left_dir.iterdir() if p.is_file() and p.name.startswith(prefix) and p.suffix == ".vtk"}
    right_files = {p.name: p for p in right_dir.iterdir() if p.is_file() and p.name.startswith(prefix) and p.suffix == ".vtk"}

    all_names = sorted(set(left_files) | set(right_files), key=lambda n: (_extract_step(n), n))
    for name in all_names:
        step = _extract_step(name)
        l = left_files.get(name)
        r = right_files.get(name)
        if l is None:
            return VtkNumericMismatch(name, kind, step, "FILE", "missing_in_left", None, None, str(r))
        if r is None:
            return VtkNumericMismatch(name, kind, step, "FILE", "missing_in_right", None, str(l), None)

        try:
            lv = _parse_vtk(l)
        except Exception as exc:
            return VtkNumericMismatch(name, kind, step, "FILE", "parse_error_left", None, str(exc), None)
        try:
            rv = _parse_vtk(r)
        except Exception as exc:
            return VtkNumericMismatch(name, kind, step, "FILE", "parse_error_right", None, None, str(exc))

        if lv.points.shape != rv.points.shape:
            return VtkNumericMismatch(name, kind, step, "POINTS", "shape_mismatch", None, str(lv.points.shape), str(rv.points.shape))
        mm = _first_allclose_mismatch(lv.points, rv.points, tol)
        if mm is not None:
            idx, d, left_val, right_val = mm
            return VtkNumericMismatch(name, kind, step, "POINTS", "value_diff", idx, left_val, right_val, d)

        if check_topology:
            if lv.cells.shape != rv.cells.shape:
                return VtkNumericMismatch(name, kind, step, "CELLS", "shape_mismatch", None, str(lv.cells.shape), str(rv.cells.shape))
            mm_cells = _first_exact_mismatch(lv.cells, rv.cells)
            if mm_cells is not None:
                idx, left_val, right_val = mm_cells
                return VtkNumericMismatch(name, kind, step, "CELLS", "value_diff", (idx,), left_val, right_val)

            if lv.cell_types.shape != rv.cell_types.shape:
                return VtkNumericMismatch(name, kind, step, "CELL_TYPES", "shape_mismatch", None, str(lv.cell_types.shape), str(rv.cell_types.shape))
            mm_ct = _first_exact_mismatch(lv.cell_types, rv.cell_types)
            if mm_ct is not None:
                idx, left_val, right_val = mm_ct
                return VtkNumericMismatch(name, kind, step, "CELL_TYPES", "value_diff", (idx,), left_val, right_val)

        left_keys = set(lv.cell_data)
        right_keys = set(rv.cell_data)
        if left_keys != right_keys:
            missing_left = sorted(right_keys - left_keys)
            missing_right = sorted(left_keys - right_keys)
            return VtkNumericMismatch(
                name,
                kind,
                step,
                "CELL_DATA",
                "scalar_name_mismatch",
                None,
                ",".join(missing_right) if missing_right else "",
                ",".join(missing_left) if missing_left else "",
            )

        for key in sorted(left_keys):
            la = lv.cell_data[key]
            ra = rv.cell_data[key]
            section = f"CELL_DATA[{key}]"
            if la.shape != ra.shape:
                return VtkNumericMismatch(name, kind, step, section, "shape_mismatch", None, str(la.shape), str(ra.shape))
            mm_cd = _first_allclose_mismatch(la, ra, tol)
            if mm_cd is not None:
                idx, d, left_val, right_val = mm_cd
                return VtkNumericMismatch(name, kind, step, section, "value_diff", idx, left_val, right_val, d)

    return None
