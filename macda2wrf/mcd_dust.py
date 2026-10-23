"""Mass-conserving MCD dust adaptation and explicit MarsWRF optical diagnostics.

MCD ``dustq`` is already kg dust / kg air; particle size does not create its
mass. Two equal-density, spherical bins reproduce mass and effective radius
only when the requested radius is between the two bin radii. Their harmonic
mass-weighted radius follows M3/M2. Optical estimates use specified constant
efficiencies, not a substitute for a spectrally resolved retrieval.

Source contracts: MCD.F90:2205-2229; Registry.EM_PLANET:515-518;
module_mp_sedim_dust.F:141-170; module_ra_mars_common.F:227-249. Radiation
receives p8w interfaces from module_radiation_driver.F:1096. Transport uses the
unweighted mass tracers, while legacy two-bin radiation weights and floors
them a second time. This module diagnoses that discrepancy without silently
changing the model or rescaling the transported mass to hide it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class DustBinParameters:
    """Sample namelist bins plus separately identified optical assumptions."""

    # 直径
    diameter1_m: float = 1.0e-6
    diameter2_m: float = 2.4e-6
    # 粒径
    qext1: float = 3.04
    qext2: float = 3.04
    # 权重
    weight1: float = 0.3
    weight2: float = 0.7
    # MCD 中的参数
    mcd_qext_visible: float = 2.4
    mcd_particle_density: float = 2500.0
    # WRF 中的参宿和
    wrf_radiation_density: float = 2650.0
    wrf_namelist_density: float = 2500.0
    # 火星常数
    gravity_m_s2: float = 3.727
    wrf_gravity_m_s2: float = 3.727
    q_floor: float = 1.0e-8
    visible_to_ir_absorption: float = 2.6

    def __post_init__(self) -> None:
        values = vars(self)
        # 值不能为无穷、负值、等于零
        for name, value in values.items():
            if not np.isfinite(value):
                raise ValueError(f"{name} must be finite")
            if name in ("weight1", "weight2", "q_floor"):
                if value < 0.0:
                    raise ValueError(f"{name} must be nonnegative")
            elif value <= 0.0:
                raise ValueError(f"{name} must be positive")
        if self.diameter1_m >= self.diameter2_m:
            raise ValueError("diameter1_m must be smaller than diameter2_m")
        if not np.isclose(self.weight1 + self.weight2, 1.0):
            raise ValueError("radiation weights must sum to one")

    @property
    def radius1_m(self) -> float:
        return 0.5 * self.diameter1_m

    @property
    def radius2_m(self) -> float:
        return 0.5 * self.diameter2_m


def _nonnegative(value: np.ndarray, name: str) -> np.ndarray:
    """Preserve NaNs as missing; reject infinities and negative finite data."""

    result = np.asarray(value, dtype=np.float64)
    if np.any(np.isinf(result)) or np.any(result < 0.0):
        raise ValueError(f"{name} must be nonnegative or NaN")
    return result


# 分箱中的有效半径
def effective_radius_from_bins(q1: np.ndarray, q2: np.ndarray, radius1_m: float, radius2_m: float) -> np.ndarray:
    """Return spherical M3/M2 radius; zero total mass has undefined radius."""

    radii = np.asarray([radius1_m, radius2_m], dtype=float)
    # 判断粒径无效值
    if np.any(~np.isfinite(radii)) or np.any(radii <= 0.0):
        raise ValueError("bin radii must be finite and positive")
    first, second = np.broadcast_arrays(_nonnegative(q1, "q1"), _nonnegative(q2, "q2"))
    total = first + second
    inverse_moment = first / radius1_m + second / radius2_m
    result = np.full(total.shape, np.nan, dtype=np.float64)
    return np.divide(total, inverse_moment, out=result, where=inverse_moment > 0.0)


def split_dust_bins(dustq: np.ndarray, reff_m: np.ndarray, parameters: DustBinParameters | None = None, strict: bool = False) -> dict[str, np.ndarray]:
    """Split physical mass and match M3/M2, with explicit range diagnostics.

    ``trc01`` and ``trc02`` are physical mass fractions, identical to ``q_bin1``
    and ``q_bin2``. Dividing these by the radiation weights would inflate the
    transported mass and is deliberately not offered as an implicit default.
    NaNs propagate. Radius clipping is recorded; strict mode rejects it.
    """

    params = parameters or DustBinParameters()
    mass, radius = np.broadcast_arrays(_nonnegative(dustq, "dustq"), np.asarray(reff_m, dtype=float))
    invalid = np.isinf(radius) | ((radius <= 0.0) & (mass > 0.0))
    if np.any(invalid):
        raise ValueError("positive dust requires a positive finite effective radius")
    clipped = ((radius < params.radius1_m) | (radius > params.radius2_m)) & (mass > 0.0)
    if strict and np.any(clipped):
        raise ValueError("MCD effective radius lies outside the configured two-bin range")
    target = np.clip(radius, params.radius1_m, params.radius2_m)
    target = np.where(mass == 0.0, params.radius1_m, target)
    denominator = 1.0 / params.radius1_m - 1.0 / params.radius2_m
    fraction = np.clip((1.0 / target - 1.0 / params.radius2_m) / denominator, 0.0, 1.0)
    first = mass * fraction
    second = mass - first
    matched = effective_radius_from_bins(first, second, params.radius1_m, params.radius2_m)
    weighted1 = np.maximum(params.q_floor, params.weight1 * first)
    weighted2 = np.maximum(params.q_floor, params.weight2 * second)
    weighted_radius = effective_radius_from_bins(weighted1, weighted2, params.radius1_m, params.radius2_m)
    return {"q_bin1": first, "q_bin2": second, "trc01": first, "trc02": second,
        "reff_requested": radius, "reff_matched": matched, "radius_clipped": clipped,
        "reff_residual": matched - radius, "radiation_q_bin1": weighted1, "radiation_q_bin2": weighted2,
        "radiation_reff": weighted_radius, "radiation_mass_residual": weighted1 + weighted2 - mass,
        "mass_residual": first + second - mass}


def pressure_layer_thickness(pressure_interfaces: np.ndarray,
    layer_count: int, axis: int = 0) -> np.ndarray:
    """Require strictly monotonic positive interfaces, in either orientation."""

    pressure = np.asarray(pressure_interfaces, dtype=np.float64)
    if pressure.ndim == 0 or not -pressure.ndim <= axis < pressure.ndim:
        raise ValueError("invalid pressure vertical axis")
    if pressure.shape[axis] != layer_count + 1:
        raise ValueError("pressure interfaces must have n_layers + 1 entries")
    if np.any(~np.isfinite(pressure)) or np.any(pressure <= 0.0):
        raise ValueError("pressure interfaces must be finite and positive")
    difference = np.diff(pressure, axis=axis)
    descending = np.all(difference < 0.0, axis=axis)
    ascending = np.all(difference > 0.0, axis=axis)
    if not np.all(descending | ascending):
        raise ValueError("every pressure column must be strictly monotonic")
    return np.abs(difference)


def column_dust_mass(dustq: np.ndarray,
    pressure_interfaces: np.ndarray, gravity_m_s2: float = 3.727, axis: int = 0) -> np.ndarray:
    """Integrate kg/kg against dp/g to kg/m2; missing layers remain missing."""

    mass = _nonnegative(dustq, "dustq")
    if mass.ndim == 0 or not -mass.ndim <= axis < mass.ndim:
        raise ValueError("dustq must have a valid vertical dimension")
    if not np.isfinite(gravity_m_s2) or gravity_m_s2 <= 0.0:
        raise ValueError("gravity_m_s2 must be finite and positive")
    thickness = pressure_layer_thickness(pressure_interfaces, mass.shape[axis], axis)
    if thickness.shape != mass.shape:
        raise ValueError("dustq and pressure-layer shapes must match exactly")
    return np.sum(mass * thickness / gravity_m_s2, axis=axis)


def column_dust_diagnostics(dustq: np.ndarray,
    reff_m: np.ndarray, pressure_interfaces: np.ndarray, parameters: DustBinParameters | None = None,
    axis: int = 0, strict: bool = False) -> dict[str, np.ndarray]:
    """Report source approximations, physical-bin optics and native WRF optics.

    The native diagnostic exactly applies the source's wt*q floors, density
    2650 and g=3.727 defaults. Its 700-Pa normalization is a historical optical
    output definition, unrelated to thermodynamic P0=610 Pa. The source visible
    estimate uses constant Q=2.4 and rho=2500 only as an explicit approximation.
    Top pressure remains positive; its omitted column is not silently filled.
    """

    params = parameters or DustBinParameters()
    bins = split_dust_bins(dustq, reff_m, params, strict)
    mass = bins["q_bin1"] + bins["q_bin2"]
    if mass.ndim == 0 or not -mass.ndim <= axis < mass.ndim:
        raise ValueError("dustq must have a valid vertical dimension")
    thickness = pressure_layer_thickness(pressure_interfaces, mass.shape[axis], axis)
    if thickness.shape != mass.shape:
        raise ValueError("dustq and pressure-layer shapes must match exactly")
    pressure = np.asarray(pressure_interfaces, dtype=float)
    surface = np.max(pressure, axis=axis)
    source_radius = np.where(mass == 0.0, params.radius1_m, bins["reff_requested"])
    source_coefficient = 0.75 * params.mcd_qext_visible / (params.mcd_particle_density * source_radius)
    source_layer_tau = source_coefficient * mass * thickness / params.gravity_m_s2
    native_coefficient = 0.75 / params.wrf_radiation_density
    unweighted = bins["q_bin1"] * params.qext1 / params.radius1_m
    unweighted += bins["q_bin2"] * params.qext2 / params.radius2_m
    weighted = bins["radiation_q_bin1"] * params.qext1 / params.radius1_m
    weighted += bins["radiation_q_bin2"] * params.qext2 / params.radius2_m
    native_layer_tau = weighted * thickness * native_coefficient / params.wrf_gravity_m_s2
    bin_layer_tau = unweighted * thickness * native_coefficient / params.wrf_gravity_m_s2
    source_tau = np.sum(source_layer_tau, axis=axis)
    bin_tau = np.sum(bin_layer_tau, axis=axis)
    native_tau = np.sum(native_layer_tau, axis=axis)
    result = dict(bins)
    result.update({"column_dust_mass": np.sum(mass * thickness / params.gravity_m_s2, axis=axis),
        "column_dust_mass_bin1": np.sum(bins["q_bin1"] * thickness / params.gravity_m_s2, axis=axis),
        "column_dust_mass_bin2": np.sum(bins["q_bin2"] * thickness / params.gravity_m_s2, axis=axis),
        "tau_visible_spherical": source_tau, "tau_visible610_spherical": source_tau * 610.0 / surface,
        "tau_ir_abs610_assumed": source_tau * 610.0 / surface / params.visible_to_ir_absorption,
        "tau_bins_visible610_unweighted": bin_tau * 610.0 / surface,
        "tau_wrf_native": native_tau, "tau_wrf_native610": native_tau * 610.0 / surface,
        "tau_wrf_native700": native_tau * 700.0 / surface,
        "tau_wrf_weighting_residual610": (native_tau - bin_tau) * 610.0 / surface})
    return result


def column_scale_dust(dustq: np.ndarray,
    background_tau: np.ndarray, target_tau: np.ndarray, axis: int = 0) -> dict[str, np.ndarray]:
    """Optional column constraint, not assimilation or vertical-profile retrieval.

    Missing targets leave background mass unchanged. Nonzero targets with zero
    background optical depth cannot determine a new profile and are rejected.
    No floor is included, so legacy WRF radiation need not close after scaling.
    """

    mass = _nonnegative(dustq, "dustq")
    if mass.ndim == 0 or not -mass.ndim <= axis < mass.ndim:
        raise ValueError("dustq must have a valid vertical dimension")
    background = _nonnegative(background_tau, "background_tau")
    target = _nonnegative(target_tau, "target_tau")
    column_shape = mass.shape[:axis] + mass.shape[axis + 1:] if axis >= 0 else np.delete(mass.shape, axis)
    column_shape = tuple(column_shape)
    background, target = np.broadcast_arrays(background, target)
    if background.shape != column_shape:
        raise ValueError("optical-depth arrays must match the nonvertical dust shape")
    valid = np.isfinite(background) & np.isfinite(target)
    if np.any(valid & (background == 0.0) & (target > 0.0)):
        raise ValueError("cannot infer nonzero dust profile from zero background optical depth")
    scale = np.ones(background.shape, dtype=float)
    np.divide(target, background, out=scale, where=valid & (background > 0.0))
    updated = mass * np.expand_dims(scale, axis=axis)
    return {"dustq": updated, "column_scale": scale, "constraint_applied": valid}
