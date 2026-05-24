# BMTF Open-Source Code (Manuscript Version)

This folder contains a cleaned, manuscript-oriented code organization for:

- Satellite data preparation and export in **Google Earth Engine (GEE)**
- **BMTF** (Bayesian Multi-source Temporal Fusion) recursive inference in Python

The goal is to make the method logic clear for paper readers and reviewers.

## Important Notes

- This package is intentionally designed for **method transparency**, not one-click reproduction.
- **No dataset or ROI assets are included**.
- You need to provide your own GEE assets, exported rasters, and sample data.
- Some scripts are simplified abstractions of the production pipeline.

## Folder Structure

```text
BMTF/
├─ LICENSE
├─ CITATION.cff
├─ .gitignore
├─ gee/
│  └─ 01_prepare_10day_dataset.js
├─ python/
│  ├─ requirements.txt
│  ├─ example_config.json
│  ├─ README.md
│  ├─ run_bmtf.py
│  └─ bmtf/
│     ├─ __init__.py
│     ├─ config.py
│     ├─ io.py
│     ├─ sar.py
│     └─ model.py
├─ examples/
│  └─ minimal_dataset/
│     ├─ README.md
│     ├─ metadata/
│     │  └─ steps_metadata.csv
│     ├─ STEP_00_YYYYMMDD_YYYYMMDD/
│     │  └─ README.md
│     └─ STEP_01_YYYYMMDD_YYYYMMDD/
│        └─ README.md
└─ CODE_MAP.md
```

## Workflow Diagram

```mermaid
flowchart LR
    A["Sentinel-2 SR/L1C + Sentinel-1 GRD (GEE)"] --> B["Quality control + TCG + S1I + 10-day windows"]
    B --> C["Export per-step GeoTIFFs and steps metadata"]
    C --> D["Local step folders (TCG/valid/nValid/S1I)"]
    D --> E["BMTF Inference: Prior + Optical + SAR (AEA)"]
    E --> F["Posterior class probabilities + MAP class"]
    F --> G["Online parameter evolution (OPE) for next step"]
```

## Recommended Workflow

1. In GEE, run [`gee/01_prepare_10day_dataset.js`](gee/01_prepare_10day_dataset.js) to export:
   - `TCG.tif`
   - `valid.tif`
   - `nValid.tif`
   - `S1I.tif`
   - `RGB.tif`
   - `steps_metadata.csv`
2. Organize local files by step folder:
   - `STEP_00_YYYYMMDD_YYYYMMDD/TCG.tif`, etc.
3. Run Python inference:
   - `python run_bmtf.py --data-root <your_data_root> --out-dir <your_output_dir> --config python/example_config.json`

## Minimal Example Directory

See [`examples/minimal_dataset`](examples/minimal_dataset/README.md) for the minimal input directory schema expected by `python/run_bmtf.py`.

## Citation Intent

This code folder corresponds to the manuscript methodology chapters:

- Section 2: Data preparation and time-window organization
- Section 3: BMTF recursive Bayesian updating + AEA + OPE

See [`CODE_MAP.md`](CODE_MAP.md) for equation-to-code mapping.
