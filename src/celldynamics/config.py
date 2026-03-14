from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ParameterFile:
    """Simple key-value parameter reader for `parameter.2dv` style files."""

    values: dict[str, str]

    @classmethod
    def load(cls, path: str | Path) -> "ParameterFile":
        p = Path(path)
        values: dict[str, str] = {}
        for line in p.read_text(encoding="utf-8").splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if "=" in s:
                k, v = s.split("=", 1)
            else:
                parts = s.split()
                if len(parts) < 2:
                    continue
                k, v = parts[0], " ".join(parts[1:])
            values[k.strip()] = v.strip()
        return cls(values=values)

    def get_float(self, key: str) -> float:
        return float(self.values[key])

    def get_int(self, key: str) -> int:
        return int(float(self.values[key]))

    def get_str(self, key: str) -> str:
        return self.values[key]
