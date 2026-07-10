"""Vertical coordinate conversion for MACDA sigma levels."""

from __future__ import annotations

import numpy as np


def sigma_pressure(psfc: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    """Return pressure on MACDA sigma levels, p(k,j,i) = psfc(j,i) * sigma(k)."""

    sigma = np.asarray(sigma, dtype=np.float64)
    psfc = np.asarray(psfc, dtype=np.float64)
    return sigma[:, None, None] * psfc[None, :, :]


def interp_sigma_to_pressure(
    field: np.ndarray,
    psfc: np.ndarray,
    sigma: np.ndarray,
    target_plev_pa: list[float] | np.ndarray,
) -> np.ndarray:
    """Interpolate a MACDA 3-D field from sigma levels to pressure levels.

    Values above the model top or below the local surface are filled by nearest
    endpoint extrapolation. This mirrors common pressure-level preprocessing
    behavior and avoids NaNs in WRF intermediate slabs.
    """

    field = np.asarray(field, dtype=np.float64)
    pressures = sigma_pressure(psfc, sigma)
    target = np.asarray(target_plev_pa, dtype=np.float64)
    out = np.empty((target.size, field.shape[-2], field.shape[-1]), dtype=np.float32)

    log_target = np.log(np.maximum(target, 1.0e-12))
    for j in range(field.shape[-2]):
        for i in range(field.shape[-1]):
            p_col = pressures[:, j, i]
            f_col = field[:, j, i]
            valid = np.isfinite(p_col) & np.isfinite(f_col) & (p_col > 0.0)
            if valid.sum() < 2:
                out[:, j, i] = np.nan
                continue
            x = np.log(p_col[valid])
            y = f_col[valid]
            order = np.argsort(x)
            x = x[order]
            y = y[order]
            out[:, j, i] = np.interp(log_target, x, y, left=y[0], right=y[-1])
    return out

