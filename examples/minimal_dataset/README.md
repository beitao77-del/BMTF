# Minimal Dataset Layout (Structure Only)

This folder demonstrates the minimal file structure expected by:

- `python/run_bmtf.py`

It is a **schema example** only. No real data is included.

## Required layout

```text
minimal_dataset/
├─ STEP_00_YYYYMMDD_YYYYMMDD/
│  ├─ TCG.tif
│  ├─ valid.tif
│  ├─ nValid.tif
│  └─ S1I.tif              # optional but recommended
├─ STEP_01_YYYYMMDD_YYYYMMDD/
│  ├─ TCG.tif
│  ├─ valid.tif
│  ├─ nValid.tif
│  └─ S1I.tif
└─ metadata/
   ├─ steps_metadata.csv
   └─ roi_out.tif          # optional
```

## Naming conventions

- Step folders should start with `STEP_`.
- Inside each step folder, expected filenames are fixed:
  - `TCG.tif`, `valid.tif`, `nValid.tif`, `S1I.tif`.
- `steps_metadata.csv` is optional but recommended for SAR timing metadata.
