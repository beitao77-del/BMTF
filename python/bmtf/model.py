from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import PipelineConfig
from .sar import OtsuResult, otsu_with_quality


def gaussian_logpdf(x: np.ndarray, mu: float, var: float, *, var_floor: float) -> np.ndarray:
    var = float(max(var, var_floor))
    return -0.5 * (((x - mu) ** 2) / var + np.log(2.0 * np.pi * var))


def softmax3(s1: np.ndarray, s2: np.ndarray, s3: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    m = np.maximum(np.maximum(s1, s2), s3)
    e1 = np.exp(s1 - m)
    e2 = np.exp(s2 - m)
    e3 = np.exp(s3 - m)
    z = e1 + e2 + e3
    return e1 / z, e2 / z, e3 / z


def _linear_clamp01(x: float, x0: float, x1: float) -> float:
    if x1 == x0:
        return 0.0
    return float(np.clip((x - x0) / (x1 - x0), 0.0, 1.0))


@dataclass
class GlobalState:
    mu1: float
    var1: float
    mu2: float
    var2: float
    mu3: float
    var3: float
    prev_p1: np.ndarray
    prev_p2: np.ndarray
    prev_p3: np.ndarray


@dataclass(frozen=True)
class SceneDiagnostics:
    qopt: float
    mode_id: int
    w_opt: float
    beta_t: float
    eta_t: float
    w_sar: float
    sar_reliable: bool
    otsu: OtsuResult | None


def markov_log_prior(
    p1: np.ndarray,
    p2: np.ndarray,
    p3: np.ndarray,
    cfg: PipelineConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    pi_t = A^T * p_{t-1}, then convert to log prior.
    """
    u0 = 1.0 / 3.0
    p1 = np.where(np.isfinite(p1), p1, u0)
    p2 = np.where(np.isfinite(p2), p2, u0)
    p3 = np.where(np.isfinite(p3), p3, u0)

    pi1 = p1 * cfg.a11 + p2 * cfg.a21 + p3 * cfg.a31
    pi2 = p1 * cfg.a12 + p2 * cfg.a22 + p3 * cfg.a32
    pi3 = p1 * cfg.a13 + p2 * cfg.a23 + p3 * cfg.a33

    pi1 = np.maximum(pi1, cfg.eps_prior)
    pi2 = np.maximum(pi2, cfg.eps_prior)
    pi3 = np.maximum(pi3, cfg.eps_prior)
    return np.log(pi1), np.log(pi2), np.log(pi3)


def _ew_update(
    prev_mu: float,
    prev_var: float,
    *,
    n: int,
    mean: float | None,
    vari: float | None,
    alpha: float,
    cfg: PipelineConfig,
) -> tuple[float, float]:
    var_floor = cfg.sigma_floor * cfg.sigma_floor
    prev_var = float(max(prev_var, var_floor))
    if mean is None or vari is None or n < cfg.min_n_update:
        return float(prev_mu), float(prev_var)

    mean = float(mean)
    vari = float(max(vari, var_floor))

    mu_cand = prev_mu * (1.0 - alpha) + mean * alpha
    shift = mean - prev_mu
    var_cand = prev_var * (1.0 - alpha) + vari * alpha + (shift * shift) * alpha * (1.0 - alpha)
    var_cand = float(max(var_cand, var_floor))

    dmu = np.clip(mu_cand - prev_mu, -cfg.mu_delta_max, cfg.mu_delta_max)
    mu_new = float(prev_mu + dmu)

    ratio = var_cand / prev_var
    ratio = float(np.clip(ratio, cfg.var_ratio_min, cfg.var_ratio_max))
    var_new = float(max(prev_var * ratio, var_floor))
    return mu_new, var_new


def infer_step(
    *,
    tcg: np.ndarray,
    valid: np.ndarray,
    nvalid: np.ndarray,
    s1i: np.ndarray | None,
    roi_mask: np.ndarray | None,
    state: GlobalState,
    cfg: PipelineConfig,
    dt_days: float | None = None,
    n_s1_eff: int | None = None,
    fallback_used: int | None = None,
) -> tuple[dict[str, np.ndarray], GlobalState, SceneDiagnostics]:
    """
    One recursive BMTF step:
    1) Prior propagation (Markov)
    2) AEA fusion (temporal + optical + SAR)
    3) Posterior update and MAP class
    4) OPE parameter update
    """
    tcg = tcg.astype(np.float32, copy=False)
    valid = (valid.astype(np.uint8, copy=False) == 1)
    nvalid = nvalid.astype(np.float32, copy=False)

    if roi_mask is None:
        in_roi = np.ones_like(valid, dtype=bool)
    else:
        in_roi = roi_mask.astype(bool, copy=False)

    opt_ok = valid & np.isfinite(tcg) & in_roi
    n_roi = int(in_roi.sum())
    qopt = float(opt_ok.sum() / n_roi) if n_roi > 0 else 0.0
    optical_absent = (qopt <= cfg.qopt_poor) or (opt_ok.sum() == 0)

    # --- SAR reliability and Otsu ---
    otsu: OtsuResult | None = None
    sar_ok = np.zeros_like(in_roi, dtype=bool)
    if s1i is not None:
        s1i = s1i.astype(np.float32, copy=False)
        sar_ok = np.isfinite(s1i) & in_roi
        if sar_ok.any():
            vals = s1i[sar_ok]
            otsu = otsu_with_quality(
                vals,
                max_buckets=cfg.otsu_max_buckets,
                min_bucket_width=cfg.hist_min_bucket_width,
                fallback_threshold=float(np.nanpercentile(vals, 90)),
                delta_min=cfg.delta_min,
            )

    dt = float(dt_days) if dt_days is not None else 0.0
    fb = int(fallback_used) if fallback_used is not None else 0

    if otsu is not None:
        q_otsu = _linear_clamp01(otsu.q, cfg.qsar_q0, cfg.qsar_q1)
        q_time = float(np.exp(-dt / max(cfg.tau_dt_days, 1e-6)))
        if n_s1_eff is not None and cfg.n_s1_ref > 0:
            q_n = float(np.clip(n_s1_eff / cfg.n_s1_ref, 0.0, 1.0))
        else:
            q_n = 1.0
        w_fallback = cfg.sar_fallback_penalty if fb == 1 else 1.0
        w_sar = float(np.clip(min(q_otsu, q_time, q_n) * w_fallback, 0.0, 1.0))
        sar_reliable = (otsu.q >= cfg.qsar_min) and (dt <= cfg.dtsar_max_days) and (w_sar >= cfg.wsar_min)
    else:
        w_sar = 0.0
        sar_reliable = False

    # --- Scene mode ---
    if not optical_absent:
        mode_id = 1 if qopt >= cfg.qopt_good else 2
    else:
        mode_id = 3 if sar_reliable else 4

    u = _linear_clamp01(qopt, cfg.qopt_poor, cfg.qopt_good)
    if mode_id == 1:
        w_opt = 1.0
        beta_t = cfg.beta_mode1
        eta_t = cfg.eta_mode1
    elif mode_id == 2:
        w_opt = 0.25 + 0.75 * u
        beta_t = cfg.beta_mode2_max - (cfg.beta_mode2_max - cfg.beta_mode2_min) * u
        eta_t = cfg.eta_mode2_max - (cfg.eta_mode2_max - cfg.eta_mode2_min) * u
    elif mode_id == 3:
        w_opt = 0.0
        beta_t = cfg.beta_poor
        eta_t = cfg.eta_mode3
    else:
        w_opt = 0.0
        beta_t = cfg.beta_poor
        eta_t = 0.0

    # --- Optical likelihood ---
    var_floor = cfg.sigma_floor * cfg.sigma_floor
    ll1_raw = gaussian_logpdf(tcg, state.mu1, state.var1, var_floor=var_floor)
    ll2_raw = gaussian_logpdf(tcg, state.mu2, state.var2, var_floor=var_floor)
    ll3_raw = gaussian_logpdf(tcg, state.mu3, state.var3, var_floor=var_floor)
    ll1_raw = np.where(np.isfinite(ll1_raw), ll1_raw, 0.0)
    ll2_raw = np.where(np.isfinite(ll2_raw), ll2_raw, 0.0)
    ll3_raw = np.where(np.isfinite(ll3_raw), ll3_raw, 0.0)

    # --- Temporal prior ---
    lp1, lp2, lp3 = markov_log_prior(state.prev_p1, state.prev_p2, state.prev_p3, cfg)

    # --- Pixel-level fallback ---
    w_opt_pix = (w_opt * opt_ok.astype(np.float32)).astype(np.float32)
    fallback_pix = (~opt_ok) & in_roi

    beta_pix = np.full_like(tcg, float(beta_t), dtype=np.float32)
    beta_pix[fallback_pix] = float(cfg.beta_poor)

    eta_pix = np.zeros_like(tcg, dtype=np.float32)
    if sar_reliable:
        eta_pix[sar_ok & opt_ok] = float(eta_t * cfg.sar_in_opt_scale)
        eta_pix[sar_ok & fallback_pix] = float(cfg.eta_mode3)

    if s1i is not None and otsu is not None:
        z = (s1i - float(otsu.threshold)) / float(otsu.delta)
        z = np.where(np.isfinite(z), z, 0.0)
        g = np.tanh(z).astype(np.float32)
    else:
        g = np.zeros_like(tcg, dtype=np.float32)

    sar_bias = g * float(w_sar) * eta_pix

    # --- AEA fused scores ---
    s1 = ll1_raw * w_opt_pix + lp1 * beta_pix - sar_bias
    s2 = ll2_raw * w_opt_pix + lp2 * beta_pix - sar_bias
    s3 = ll3_raw * w_opt_pix + lp3 * beta_pix + sar_bias

    p1_upd, p2_upd, p3_upd = softmax3(s1, s2, s3)
    cls = (np.argmax(np.stack([p1_upd, p2_upd, p3_upd], axis=0), axis=0) + 1).astype(np.int16)

    # Output confidence support: optical=1, SAR-only=w_sar, none=0.
    r_out = np.where(opt_ok, 1.0, np.where(sar_ok, float(w_sar), 0.0)).astype(np.float32)
    u0 = 1.0 / 3.0
    p1 = ((1.0 - r_out) * u0 + r_out * p1_upd).astype(np.float32)
    p2 = ((1.0 - r_out) * u0 + r_out * p2_upd).astype(np.float32)
    p3 = ((1.0 - r_out) * u0 + r_out * p3_upd).astype(np.float32)
    top1 = np.maximum(np.maximum(p1, p2), p3).astype(np.float32)

    cls = np.where(in_roi, cls, 0).astype(np.int16)
    p1 = np.where(in_roi, p1, np.nan).astype(np.float32)
    p2 = np.where(in_roi, p2, np.nan).astype(np.float32)
    p3 = np.where(in_roi, p3, np.nan).astype(np.float32)
    top1 = np.where(in_roi, top1, np.nan).astype(np.float32)

    # --- OPE: high-confidence update samples ---
    allow_update_scene = (mode_id == 1) or (mode_id == 2 and cfg.allow_update_mode2)
    update_mask = (top1 >= cfg.t_update) & opt_ok & (nvalid >= cfg.nvalid_upd) & allow_update_scene

    def _mean_var_count(mask: np.ndarray) -> tuple[float | None, float | None, int]:
        idx = np.where(mask)
        if idx[0].size == 0:
            return None, None, 0
        v = tcg[idx].astype(np.float64, copy=False)
        return float(v.mean()), float(v.var(ddof=0)), int(v.size)

    m1, v1, n1 = _mean_var_count(update_mask & (cls == 1))
    m2, v2, n2 = _mean_var_count(update_mask & (cls == 2))
    m3, v3, n3 = _mean_var_count(update_mask & (cls == 3))

    alpha = float(
        np.clip(
            cfg.alpha_min + (cfg.alpha_max - cfg.alpha_min) * np.clip(qopt, 0.0, 1.0),
            cfg.alpha_min,
            cfg.alpha_max,
        )
    )

    mu1n, var1n = _ew_update(state.mu1, state.var1, n=n1, mean=m1, vari=v1, alpha=alpha, cfg=cfg)
    mu2n, var2n = _ew_update(state.mu2, state.var2, n=n2, mean=m2, vari=v2, alpha=alpha, cfg=cfg)
    mu3n, var3n = _ew_update(state.mu3, state.var3, n=n3, mean=m3, vari=v3, alpha=alpha, cfg=cfg)

    # Hidden state propagation can be conservative in optical gaps.
    r_state = np.where(opt_ok, 1.0, np.where(sar_ok, float(w_sar) * float(cfg.state_sar_scale), 0.0)).astype(np.float32)
    p1_state = ((1.0 - r_state) * state.prev_p1 + r_state * p1_upd).astype(np.float32)
    p2_state = ((1.0 - r_state) * state.prev_p2 + r_state * p2_upd).astype(np.float32)
    p3_state = ((1.0 - r_state) * state.prev_p3 + r_state * p3_upd).astype(np.float32)
    p1_state = np.where(in_roi, p1_state, np.nan).astype(np.float32)
    p2_state = np.where(in_roi, p2_state, np.nan).astype(np.float32)
    p3_state = np.where(in_roi, p3_state, np.nan).astype(np.float32)

    next_state = GlobalState(
        mu1=mu1n,
        var1=var1n,
        mu2=mu2n,
        var2=var2n,
        mu3=mu3n,
        var3=var3n,
        prev_p1=p1_state,
        prev_p2=p2_state,
        prev_p3=p3_state,
    )

    outputs = {
        "class": cls,
        "top1_post": top1,
        "p1": p1,
        "p2": p2,
        "p3": p3,
        "update_mask": update_mask.astype(np.uint8),
    }

    diag = SceneDiagnostics(
        qopt=float(qopt),
        mode_id=int(mode_id),
        w_opt=float(w_opt),
        beta_t=float(beta_t),
        eta_t=float(eta_t),
        w_sar=float(w_sar),
        sar_reliable=bool(sar_reliable),
        otsu=otsu,
    )
    return outputs, next_state, diag

