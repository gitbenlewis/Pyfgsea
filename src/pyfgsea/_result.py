"""Small plotting result wrapper."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class PlotResult:
    """Container returned by plotting helpers."""

    figure: Any
    axes: Any
    data: Any = None

    def save(self, path: str | Path, **kwargs: Any) -> None:
        self.figure.savefig(path, **kwargs)
