"""Regular latitude-longitude grid helpers."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import RegularGridInterpolator


@dataclass(frozen=True)
class RegularLatLonGrid:
    lat_start: float
    lat_end: float
    nlat: int
    lon_start: float
    lon_end: float
    nlon: int
    radius_km: float

    @property
    def lats(self) -> np.ndarray:
        return np.linspace(self.lat_start, self.lat_end, self.nlat)

    @property
    def lons(self) -> np.ndarray:
        return np.linspace(self.lon_start, self.lon_end, self.nlon)

    @property
    def deltlat(self) -> float:
        return (self.lat_end - self.lat_start) / (self.nlat - 1)

    @property
    def deltlon(self) -> float:
        return (self.lon_end - self.lon_start) / (self.nlon - 1)


def orient_lat_lon(
    arr: np.ndarray,
    src_lats: np.ndarray,
    src_lons: np.ndarray,
    lon_convention: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return data with ascending latitude and monotonic target longitudes."""

    data = np.asarray(arr)
    lats = np.asarray(src_lats, dtype=np.float64)
    lons = np.asarray(src_lons, dtype=np.float64)

    if lats[0] > lats[-1]:
        lats = lats[::-1]
        data = np.flip(data, axis=-2)

    if lon_convention == "0_360":
        lons = np.mod(lons, 360.0)
    elif lon_convention in {"native", "-180_180"}:
        pass
    else:
        raise ValueError(f"Unsupported lon_convention: {lon_convention}")

    order = np.argsort(lons)
    lons = lons[order]
    data = np.take(data, order, axis=-1)
    return data, lats, lons


def horizontal_interp_2d(
    slab: np.ndarray,
    src_lats: np.ndarray,
    src_lons: np.ndarray,
    target_grid: RegularLatLonGrid,
) -> np.ndarray:
    """Interpolate one 2-D slab onto the target regular lat-lon grid."""

    src_lats = np.asarray(src_lats, dtype=np.float64)
    src_lons = np.asarray(src_lons, dtype=np.float64)
    slab = np.asarray(slab, dtype=np.float64)

    if (
        slab.shape == (target_grid.nlat, target_grid.nlon)
        and np.allclose(src_lats, target_grid.lats)
        and np.allclose(src_lons, target_grid.lons)
    ):
        return slab.astype(np.float32, copy=False)

    interp = RegularGridInterpolator(
        (src_lats, src_lons),
        slab,
        bounds_error=False,
        fill_value=None,
    )
    lat2d, lon2d = np.meshgrid(target_grid.lats, target_grid.lons, indexing="ij")
    points = np.column_stack([lat2d.ravel(), lon2d.ravel()])
    out = interp(points).reshape(target_grid.nlat, target_grid.nlon)
    return out.astype(np.float32)

