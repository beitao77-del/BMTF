from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class OtsuResult:
    threshold: float
    q: float
    sigma_t2: float
    sigma_b2: float
    p90: float
    delta: float


def otsu_with_quality(
    values: np.ndarray,
    *,
    max_buckets: int,
    min_bucket_width: float,
    fallback_threshold: float,
    delta_min: float,
) -> OtsuResult:
    """
    Otsu threshold on 1D SAR indicator values.
    q = sigma_b^2 / sigma_t^2, used as SAR separability quality.
    """
    v = values.astype(np.float64, copy=False)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return OtsuResult(
            threshold=float(fallback_threshold),
            q=0.0,
            sigma_t2=0.0,
            sigma_b2=0.0,
            p90=float("nan"),
            delta=float(delta_min),
        )

    p90 = float(np.percentile(v, 90))
    vmin = float(np.min(v))
    vmax = float(np.max(v))
    if (not np.isfinite(vmin)) or (not np.isfinite(vmax)) or vmax <= vmin:
        thr = float(fallback_threshold)
        return OtsuResult(threshold=thr, q=0.0, sigma_t2=0.0, sigma_b2=0.0, p90=p90, delta=max(p90 - thr, delta_min))

    max_bins_by_width = int(np.ceil((vmax - vmin) / max(min_bucket_width, 1e-6)))
    bins = int(max(8, min(max_buckets, max_bins_by_width)))
    hist, bin_edges = np.histogram(v, bins=bins, range=(vmin, vmax))
    if hist.sum() == 0:
        thr = float(fallback_threshold)
        return OtsuResult(threshold=thr, q=0.0, sigma_t2=0.0, sigma_b2=0.0, p90=p90, delta=max(p90 - thr, delta_min))

    centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    w = hist.astype(np.float64)
    w = w / w.sum()

    mu_t = float(np.sum(w * centers))
    sigma_t2 = float(np.sum(w * (centers - mu_t) ** 2))
    if sigma_t2 <= 0:
        thr = float(fallback_threshold)
        return OtsuResult(threshold=thr, q=0.0, sigma_t2=sigma_t2, sigma_b2=0.0, p90=p90, delta=max(p90 - thr, delta_min))

    omega = np.cumsum(w)
    mu = np.cumsum(w * centers)
    denom = omega * (1.0 - omega)
    valid = denom > 1e-12
    sigma_b2 = np.zeros_like(centers)
    sigma_b2[valid] = (mu_t * omega[valid] - mu[valid]) ** 2 / denom[valid]

    k = int(np.argmax(sigma_b2))
    thr = float(centers[k])
    sigma_b2_max = float(sigma_b2[k])
    q = float(np.clip(sigma_b2_max / sigma_t2, 0.0, 1.0))
    delta = max(p90 - thr, float(delta_min))

    return OtsuResult(threshold=thr, q=q, sigma_t2=sigma_t2, sigma_b2=sigma_b2_max, p90=p90, delta=delta)

