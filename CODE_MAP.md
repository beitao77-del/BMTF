# Code-to-Method Map

This file links manuscript sections (Chapter 2 and 3) to the code in this folder.

## Section 2.2 Satellite Data (GEE preprocessing/export)

- Script: [01_prepare_10day_dataset.js](gee/01_prepare_10day_dataset.js)
- Key logic:
  - Sentinel-2 quality masking (`qualityMaskL2A`)
  - TCG computation from S2 L1C (`tcgFromL1C`)
  - Sentinel-1 linear-domain Lee filtering (`leeFilterLinear`)
  - `S1I = VV_dB + VH_dB`
  - 10-day window aggregation and export

## Section 3.2 Recursive Bayesian Updating

- File: [model.py](python/bmtf/model.py)
- Functions:
  - `gaussian_logpdf` (optical likelihood)
  - `markov_log_prior` (temporal prior propagation)
  - `infer_step` (posterior update and MAP classification)

## Section 3.3 Adaptive Evidence Aggregation (AEA)

- File: [model.py](python/bmtf/model.py)
- Components inside `infer_step`:
  - Scene-level mode selection (Mode1/2/3/4)
  - Adaptive weighting of optical/prior/SAR evidence
  - Signed SAR bias integration
- File: [sar.py](python/bmtf/sar.py)
  - Otsu threshold + quality metric for SAR reliability

## Section 3.4 Online Parameter Evolution (OPE)

- File: [model.py](python/bmtf/model.py)
- Functions:
  - `_ew_update` (constrained exponential-forgetting update)
  - Reliable high-confidence sample gating inside `infer_step`

## Section 3.5 Accuracy Assessment (implementation entry point)

- File: [run_bmtf.py](python/run_bmtf.py)
- Provides:
  - Sequential inference over 10-day windows
  - Output map writing (`class`, `top1_post`, `p1/p2/p3`)
  - Lightweight run log (`run_steps.csv`)
