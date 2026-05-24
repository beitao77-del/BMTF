# Python Inference Module

This folder contains a compact BMTF implementation focused on method clarity.

## Expected Input Layout

```text
<data_root>/
├─ STEP_00_YYYYMMDD_YYYYMMDD/
│  ├─ TCG.tif
│  ├─ valid.tif
│  ├─ nValid.tif
│  ├─ S1I.tif           # optional but recommended
│  └─ RGB.tif           # optional
├─ STEP_01_...
└─ metadata/
   ├─ steps_metadata.csv
   └─ roi_out.tif       # optional
```

## Run

```bash
python run_bmtf.py \
  --data-root /path/to/data_root \
  --out-dir /path/to/output \
  --config example_config.json
```

Outputs include per-step `class.tif`, `top1_post.tif`, `p1.tif`, `p2.tif`, `p3.tif`, plus `_analysis/run_steps.csv`.

For a minimal directory schema example, see:

- `../examples/minimal_dataset/`
