from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Tolerance:
    rtol: float = 1e-6
    atol: float = 1e-8


def compare_arrays(
    actual: np.ndarray,
    expected: np.ndarray,
    tol: Tolerance = Tolerance(),
    label: str = "array",
) -> None:
    if actual.shape != expected.shape:
        raise AssertionError(f"{label}: shape mismatch {actual.shape} != {expected.shape}")

    if not np.allclose(actual, expected, rtol=tol.rtol, atol=tol.atol):
        diff = np.abs(actual - expected)
        i = np.unravel_index(np.argmax(diff), diff.shape)
        raise AssertionError(
            f"{label}: max diff {diff[i]:.3e} at index {i}; "
            f"actual={actual[i]!r}, expected={expected[i]!r}, "
            f"rtol={tol.rtol}, atol={tol.atol}"
        )
