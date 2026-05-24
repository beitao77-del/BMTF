from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

try:
    import rasterio
except Exception as exc:  # pragma: no cover
    rasterio = None  # type: ignore[assignment]
    _RASTERIO_IMPORT_ERROR = exc


@dataclass(frozen=True)
class StepPaths:
    step_dir: Path
    tcg: Path
    valid: Path
    nvalid: Path
    s1i: Path
    rgb: Path


@dataclass(frozen=True)
class StepMeta:
    idx: int
    label: str
    s2_center_ms: int | None = None
    s1_center_ms: int | None = None
    n_s1_eff: int | None = None
    fallback_used: int | None = None


def _require_rasterio() -> Any:
    if rasterio is None:  # pragma: no cover
        raise RuntimeError("rasterio is required for GeoTIFF I/O.") from _RASTERIO_IMPORT_ERROR
    return rasterio


def list_step_dirs(data_root: Path) -> list[Path]:
    step_dirs = [p for p in data_root.iterdir() if p.is_dir() and p.name.lower().startswith("step_")]

    def _key(p: Path) -> tuple[int, str]:
        m = re.match(r"step_(\d+)_", p.name, flags=re.IGNORECASE)
        return (int(m.group(1)), p.name) if m else (10**9, p.name)

    return sorted(step_dirs, key=_key)


def resolve_step_paths(step_dir: Path) -> StepPaths:
    return StepPaths(
        step_dir=step_dir,
        tcg=step_dir / "TCG.tif",
        valid=step_dir / "valid.tif",
        nvalid=step_dir / "nValid.tif",
        s1i=step_dir / "S1I.tif",
        rgb=step_dir / "RGB.tif",
    )


def load_steps_csv(data_root: Path) -> dict[int, StepMeta]:
    """
    Expected path: data_root/metadata/steps_metadata.csv
    """
    csv_path = data_root / "metadata" / "steps_metadata.csv"
    if not csv_path.exists():
        return {}

    def _opt_int(row: dict[str, str], key: str) -> int | None:
        v = row.get(key, "")
        if v == "":
            return None
        try:
            return int(float(v))
        except Exception:
            return None

    out: dict[int, StepMeta] = {}
    with csv_path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            idx = int(row["idx"])
            out[idx] = StepMeta(
                idx=idx,
                label=row.get("label", f"STEP_{idx:02d}"),
                s2_center_ms=_opt_int(row, "s2_center_ms"),
                s1_center_ms=_opt_int(row, "s1_center_ms"),
                n_s1_eff=_opt_int(row, "nS1"),
                fallback_used=_opt_int(row, "fallback_used"),
            )
    return out


def read_raster(path: Path, masked: bool = True) -> tuple[np.ndarray, dict[str, Any]]:
    rio = _require_rasterio()
    with rio.open(path) as ds:
        arr = ds.read(1, masked=masked)
        profile = ds.profile.copy()
    return arr, profile


def write_raster(path: Path, arr: np.ndarray, profile: dict[str, Any]) -> None:
    rio = _require_rasterio()
    path.parent.mkdir(parents=True, exist_ok=True)
    out_profile = profile.copy()

    if arr.ndim == 2:
        out_profile.update(count=1)
    elif arr.ndim == 3:
        out_profile.update(count=arr.shape[0])
    else:
        raise ValueError(f"Unsupported array shape: {arr.shape}")

    out_profile.pop("blockxsize", None)
    out_profile.pop("blockysize", None)
    out_profile.pop("tiled", None)

    with rio.open(path, "w", **out_profile) as ds:
        if arr.ndim == 2:
            ds.write(arr, 1)
        else:
            ds.write(arr)


def read_roi_mask(data_root: Path) -> np.ndarray | None:
    """
    Optional ROI mask path:
    data_root/metadata/roi_out.tif
    """
    roi_path = data_root / "metadata" / "roi_out.tif"
    if not roi_path.exists():
        return None
    roi, _ = read_raster(roi_path, masked=False)
    return (roi > 0).astype(np.uint8)

