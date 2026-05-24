from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np

from bmtf.config import PipelineConfig
from bmtf.io import list_step_dirs, load_steps_csv, read_raster, read_roi_mask, resolve_step_paths, write_raster
from bmtf.model import GlobalState, infer_step


def _init_state(shape: tuple[int, int], cfg: PipelineConfig) -> GlobalState:
    p = np.full(shape, 1.0 / 3.0, dtype=np.float32)
    return GlobalState(
        mu1=float(cfg.mu1_init),
        var1=float(cfg.var1_init),
        mu2=float(cfg.mu2_init),
        var2=float(cfg.var2_init),
        mu3=float(cfg.mu3_init),
        var3=float(cfg.var3_init),
        prev_p1=p.copy(),
        prev_p2=p.copy(),
        prev_p3=p.copy(),
    )


def _step_index(step_name: str) -> int | None:
    parts = step_name.split("_")
    if len(parts) < 2:
        return None
    try:
        return int(parts[1])
    except Exception:
        return None


def run_sequence(
    *,
    data_root: Path,
    out_dir: Path,
    cfg: PipelineConfig,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    analysis_dir = out_dir / "_analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    (analysis_dir / "run_config.json").write_text(json.dumps(asdict(cfg), ensure_ascii=False, indent=2), encoding="utf-8")

    step_dirs = list_step_dirs(data_root)
    if not step_dirs:
        raise RuntimeError(f"No STEP_* folders found under: {data_root}")

    steps_meta = load_steps_csv(data_root)
    first_paths = resolve_step_paths(step_dirs[0])
    tcg0, profile = read_raster(first_paths.tcg, masked=False)
    state = _init_state(tcg0.shape, cfg)
    roi_mask = read_roi_mask(data_root)

    run_csv = analysis_dir / "run_steps.csv"
    with run_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "step",
                "mode_id",
                "qopt",
                "w_opt",
                "beta_t",
                "eta_t",
                "w_sar",
                "sar_reliable",
                "dt_days",
                "n_s1_eff",
                "fallback_used",
                "mu1",
                "var1",
                "mu2",
                "var2",
                "mu3",
                "var3",
            ]
        )

        for step_dir in step_dirs:
            paths = resolve_step_paths(step_dir)
            tcg, profile = read_raster(paths.tcg, masked=False)
            valid, _ = read_raster(paths.valid, masked=False)
            nvalid, _ = read_raster(paths.nvalid, masked=False)
            s1i = None
            if paths.s1i.exists():
                s1i, _ = read_raster(paths.s1i, masked=False)

            idx = _step_index(step_dir.name)
            dt_days = None
            n_s1_eff = None
            fallback_used = None
            if idx is not None and idx in steps_meta:
                m = steps_meta[idx]
                if m.s1_center_ms is not None and m.s2_center_ms is not None:
                    dt_days = abs(m.s1_center_ms - m.s2_center_ms) / (1000 * 60 * 60 * 24)
                n_s1_eff = m.n_s1_eff
                fallback_used = m.fallback_used

            outputs, state, diag = infer_step(
                tcg=tcg,
                valid=valid,
                nvalid=nvalid,
                s1i=s1i,
                roi_mask=roi_mask,
                state=state,
                cfg=cfg,
                dt_days=dt_days,
                n_s1_eff=n_s1_eff,
                fallback_used=fallback_used,
            )

            step_out = out_dir / step_dir.name
            step_out.mkdir(parents=True, exist_ok=True)

            class_profile = profile.copy()
            class_profile.update(dtype="int16", nodata=0)
            write_raster(step_out / "class.tif", outputs["class"].astype(np.int16), class_profile)

            float_profile = profile.copy()
            float_profile.update(dtype="float32", nodata=np.nan)
            write_raster(step_out / "top1_post.tif", outputs["top1_post"].astype(np.float32), float_profile)
            write_raster(step_out / "p1.tif", outputs["p1"].astype(np.float32), float_profile)
            write_raster(step_out / "p2.tif", outputs["p2"].astype(np.float32), float_profile)
            write_raster(step_out / "p3.tif", outputs["p3"].astype(np.float32), float_profile)

            w.writerow(
                [
                    step_dir.name,
                    diag.mode_id,
                    f"{diag.qopt:.6f}",
                    f"{diag.w_opt:.6f}",
                    f"{diag.beta_t:.6f}",
                    f"{diag.eta_t:.6f}",
                    f"{diag.w_sar:.6f}",
                    int(diag.sar_reliable),
                    (f"{dt_days:.6f}" if dt_days is not None else ""),
                    (int(n_s1_eff) if n_s1_eff is not None else ""),
                    (int(fallback_used) if fallback_used is not None else ""),
                    f"{state.mu1:.8f}",
                    f"{state.var1:.10f}",
                    f"{state.mu2:.8f}",
                    f"{state.var2:.10f}",
                    f"{state.mu3:.8f}",
                    f"{state.var3:.10f}",
                ]
            )

            print(
                f"{step_dir.name}: mode={diag.mode_id} qopt={diag.qopt:.3f} "
                f"wOpt={diag.w_opt:.3f} beta={diag.beta_t:.3f} eta={diag.eta_t:.3f} "
                f"wSar={diag.w_sar:.3f} mu=[{state.mu1:.4f},{state.mu2:.4f},{state.mu3:.4f}]"
            )

    print(f"Wrote: {run_csv}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="run_bmtf")
    ap.add_argument("--data-root", required=True, help="Root folder containing STEP_* and metadata/")
    ap.add_argument("--out-dir", required=True, help="Output folder for BMTF results")
    ap.add_argument("--config", default=None, help="Optional JSON config to override defaults")
    args = ap.parse_args(argv)

    cfg = PipelineConfig()
    if args.config:
        with Path(args.config).open("r", encoding="utf-8") as f:
            cfg_dict = json.load(f)
        allowed = set(cfg.__dataclass_fields__.keys())
        cfg = replace(cfg, **{k: v for k, v in cfg_dict.items() if k in allowed})

    run_sequence(data_root=Path(args.data_root), out_dir=Path(args.out_dir), cfg=cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

