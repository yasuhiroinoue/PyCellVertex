from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
import re


@dataclass(frozen=True)
class FileMismatch:
    filename: str
    reason: str
    left: str | None
    right: str | None


_STEP_RE = re.compile(r"(\d+)")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _extract_step(name: str) -> int:
    m = _STEP_RE.findall(name)
    if not m:
        return -1
    return int(m[-1])


def first_mismatch(left_dir: Path, right_dir: Path) -> FileMismatch | None:
    left_files = {p.name: p for p in left_dir.iterdir() if p.is_file()}
    right_files = {p.name: p for p in right_dir.iterdir() if p.is_file()}

    all_names = sorted(set(left_files) | set(right_files))
    for name in all_names:
        l = left_files.get(name)
        r = right_files.get(name)
        if l is None:
            return FileMismatch(name, "missing_in_left", None, str(r))
        if r is None:
            return FileMismatch(name, "missing_in_right", str(l), None)

        hl = sha256_file(l)
        hr = sha256_file(r)
        if hl != hr:
            return FileMismatch(name, "content_diff", hl, hr)

    return None


def first_mismatch_by_step(left_dir: Path, right_dir: Path, suffix: str = ".vtk") -> FileMismatch | None:
    """Find first mismatch ordered by inferred timestep from filename."""
    left_files = {p.name: p for p in left_dir.iterdir() if p.is_file() and p.name.endswith(suffix)}
    right_files = {p.name: p for p in right_dir.iterdir() if p.is_file() and p.name.endswith(suffix)}

    all_names = sorted(set(left_files) | set(right_files), key=lambda n: (_extract_step(n), n))
    for name in all_names:
        l = left_files.get(name)
        r = right_files.get(name)
        if l is None:
            return FileMismatch(name, "missing_in_left", None, str(r))
        if r is None:
            return FileMismatch(name, "missing_in_right", str(l), None)

        hl = sha256_file(l)
        hr = sha256_file(r)
        if hl != hr:
            return FileMismatch(name, "content_diff", hl, hr)

    return None
