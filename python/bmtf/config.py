from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PipelineConfig:
    # Data range
    tcg_min: float = -0.12
    tcg_max: float = 0.35

    # Initial optical likelihood parameters: N(mu_c, var_c)
    mu1_init: float = -0.02
    var1_init: float = 0.0005
    mu2_init: float = -0.05
    var2_init: float = 0.00008
    mu3_init: float = 0.06
    var3_init: float = 0.0013

    # OPE update gate
    t_update: float = 0.80
    min_n_update: int = 5000
    nvalid_upd: int = 3
    allow_update_mode2: bool = False

    # Numerical stability
    sigma_floor: float = 1e-3
    eps_prior: float = 1e-6

    # OPE hard constraints
    mu_delta_max: float = 0.010
    var_ratio_min: float = 0.80
    var_ratio_max: float = 1.20

    # AEA mode thresholds
    qopt_good: float = 0.90
    qopt_poor: float = 0.60

    # Prior weights (beta)
    beta_mode1: float = 1.0
    beta_mode2_min: float = 1.1
    beta_mode2_max: float = 1.8
    beta_poor: float = 2.0

    # SAR reliability
    qsar_q0: float = 0.15
    qsar_q1: float = 0.45
    tau_dt_days: float = 7.0
    n_s1_ref: float = 3.0
    sar_fallback_penalty: float = 0.55
    wsar_min: float = 0.15
    qsar_min: float = 0.15
    dtsar_max_days: float = 12.0

    # SAR evidence strength (eta)
    eta_mode1: float = 1.0
    eta_mode2_min: float = 1.2
    eta_mode2_max: float = 1.8
    eta_mode3: float = 2.0

    # How much SAR enters optical-valid pixels vs optical-gap pixels
    sar_in_opt_scale: float = 0.2
    state_sar_scale: float = 0.0

    # Otsu histogram settings
    otsu_max_buckets: int = 128
    hist_min_bucket_width: float = 0.02
    delta_min: float = 0.01

    # Markov transition matrix A
    a11: float = 0.90
    a12: float = 0.06
    a13: float = 0.04
    a21: float = 0.05
    a22: float = 0.92
    a23: float = 0.03
    a31: float = 0.03
    a32: float = 0.04
    a33: float = 0.93

    # OPE scene-level alpha schedule
    alpha_min: float = 0.08
    alpha_max: float = 0.32

