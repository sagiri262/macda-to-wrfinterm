"""MCD hybrid pressure, hydrostatic height and conservative vertical mapping.

All vertical arrays are ordered surface to top and have shape (level, lat, lon).
Pressure is Pa. Orography is geometric height above the local areoid in metres;
the optional ``areoid_m`` is an areoid RADIUS, not an elevation offset.

This module does not generate a WRF initial condition. WPS receives pressure and
geopotential height; real.exe constructs the model's actual dry hydrostatic eta
coordinate. ``pressure_to_eta`` is only a total-pressure diagnostic analogue.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


MARS_GRAVITY_M_S2 = 3.727
MARS_DRY_GAS_CONSTANT = 188.92
MARS_REFERENCE_PRESSURE_PA = 610.0
MARS_REFERENCE_RADIUS_M = 3389920.0


@dataclass(frozen=True)
class InterfaceResult:
    """Approximate interfaces and an auditable account of their reconstruction."""

    pressure_pa: np.ndarray
    method: str
    top_boundary: str


@dataclass(frozen=True)
class HeightResult:
    """Separate geopotential, geometric and radius-reference heights in metres."""

    geopotential_height_areoid_m: np.ndarray
    geometric_height_areoid_m: np.ndarray
    geometric_height_reference_m: np.ndarray | None
    interface_geopotential_height_areoid_m: np.ndarray
    interface_geometric_height_areoid_m: np.ndarray
    effective_gas_constant: np.ndarray
    method: str


def _finite_array(values: np.ndarray, name: str,
    dimensions: int | None = None) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if dimensions is not None and array.ndim != dimensions:
        raise ValueError(f"{name} must have {dimensions} dimensions")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains missing or nonfinite values")
    return array


def _positive_array(values: np.ndarray, name: str,
    dimensions: int | None = None) -> np.ndarray:
    array = _finite_array(values, name, dimensions)
    if np.any(array <= 0.0):
        raise ValueError(f"{name} must be strictly positive")
    return array


def _descending_pressure(values: np.ndarray,
    name: str) -> np.ndarray:
    array = _positive_array(values, name, 3)
    if array.shape[0] < 2 or np.any(np.diff(array, axis=0) >= 0.0):
        raise ValueError(f"{name} must be strictly decreasing surface to top")
    return array


def hybrid_pressure(surface_pressure_pa: np.ndarray,
    aps_pa: np.ndarray, bps: np.ndarray) -> np.ndarray:
    """Return MCD midlayer pressure A(k) + B(k) * ps(lat, lon)."""
    surface = _positive_array(surface_pressure_pa, "surface pressure", 2)
    pressure_term = _finite_array(aps_pa, "aps", 1)
    sigma_term = _finite_array(bps, "bps", 1)
    if pressure_term.shape != sigma_term.shape:
        raise ValueError("aps and bps must have equal lengths")
    pressure = pressure_term[:, None, None] + sigma_term[:, None, None] * surface[None, :, :]
    pressure = _descending_pressure(pressure, "hybrid pressure")
    if np.any(pressure[0] >= surface):
        raise ValueError("lowest pressure centre must lie strictly above the surface")
    return pressure


def pressure_interfaces(pressure_center_pa: np.ndarray,
    surface_pressure_pa: np.ndarray, top_pressure_pa: float | None = None) -> InterfaceResult:
    """Infer approximate interfaces from centres, never claim native PCM bounds.

    Interior boundaries are log-pressure midpoints. The bottom boundary equals
    ps. Without an explicit top, reflect the final centre about its lower
    interface in log pressure; a zero-pressure top cannot have a finite height.
    These bounds conserve a chosen discrete column, not the unknown native PCM
    finite-volume masses. The inference must be retained in output metadata.
    """
    centres = _descending_pressure(pressure_center_pa, "pressure centres")
    surface = _positive_array(surface_pressure_pa, "surface pressure", 2)
    if surface.shape != centres.shape[1:]:
        raise ValueError("surface pressure and pressure centres have incompatible horizontal shapes")
    if np.any(centres[0] >= surface):
        raise ValueError("lowest pressure centre must lie strictly above the surface")
    interfaces = np.empty((centres.shape[0] + 1,) + surface.shape, dtype=np.float64)
    interfaces[0] = surface
    interfaces[1:-1] = np.exp(0.5 * (np.log(centres[:-1]) + np.log(centres[1:])))
    if top_pressure_pa is None:
        interfaces[-1] = np.exp(2.0 * np.log(centres[-1]) - np.log(interfaces[-2]))
        top_method = "approximate log-symmetric top extrapolation"
    else:
        top = float(top_pressure_pa)
        if not np.isfinite(top) or top <= 0.0 or np.any(top >= centres[-1]):
            raise ValueError("explicit top pressure must be positive and above all top centres")
        interfaces[-1] = top
        top_method = "explicit top pressure"
    _descending_pressure(interfaces, "pressure interfaces")
    return InterfaceResult(interfaces, "approximate geometric midpoints; bottom=ps", top_method)


def layer_pressure_thickness(pressure_interface_pa: np.ndarray) -> np.ndarray:
    """Return positive pressure thickness, including explicit top/bottom bounds."""
    interfaces = _descending_pressure(pressure_interface_pa, "pressure interfaces")
    return interfaces[:-1] - interfaces[1:]


def hydrostatic_height(pressure_center_pa: np.ndarray,
    pressure_interface_pa: np.ndarray, temperature_k: np.ndarray,
    surface_orography_m: np.ndarray, areoid_m: np.ndarray | None = None,
    density_kg_m3: np.ndarray | None = None, gas_constant: float = MARS_DRY_GAS_CONSTANT,
    gravity: float = MARS_GRAVITY_M_S2, geometry: str = "constant_gravity",
    reference_radius_m: float = MARS_REFERENCE_RADIUS_M) -> HeightResult:
    """Integrate hydrostatic geopotential with piecewise-isothermal layers.

    R=p/(rho*T) follows the MCD accessor when density is available. Otherwise
    a declared constant Mars Rd is used. This finite-layer integration is an
    explicit approximation, not a byte-for-byte reproduction of MCD's accessor.

    In constant_gravity mode geometric and geopotential heights coincide. The
    optional inverse_square mode assumes spherical gravity g(z)=g0*r0^2/(r0+z)^2
    about the local areoid radius r0; it maps geopotential height H to geometric
    z=H*r0/(r0-H). Terrain is never moved and ps is never adjusted by this routine.
    """
    centres = _descending_pressure(pressure_center_pa, "pressure centres")
    interfaces = _descending_pressure(pressure_interface_pa, "pressure interfaces")
    temperature = _positive_array(temperature_k, "temperature", 3)
    orography = _finite_array(surface_orography_m, "orography", 2)
    if temperature.shape != centres.shape or interfaces.shape != (centres.shape[0] + 1,) + centres.shape[1:]:
        raise ValueError("temperature, pressure centres and interfaces must describe the same layers")
    if orography.shape != centres.shape[1:]:
        raise ValueError("orography and pressure have incompatible horizontal shapes")
    if np.any(centres >= interfaces[:-1]) or np.any(centres <= interfaces[1:]):
        raise ValueError("each pressure centre must lie strictly inside its pressure interfaces")
    if not np.isfinite(gravity) or gravity <= 0.0:
        raise ValueError("gravity must be positive and finite")
    if not np.isfinite(gas_constant) or gas_constant <= 0.0:
        raise ValueError("gas constant must be positive and finite")
    if not np.isfinite(reference_radius_m) or reference_radius_m <= 0.0:
        raise ValueError("reference radius must be positive and finite")
    if geometry not in ("constant_gravity", "inverse_square"):
        raise ValueError("unknown gravity geometry")
    if density_kg_m3 is None:
        effective_r = np.full_like(centres, gas_constant)
        gas_method = f"constant R={gas_constant:g} J kg-1 K-1"
    else:
        density = _positive_array(density_kg_m3, "density", 3)
        if density.shape != centres.shape:
            raise ValueError("density and pressure must have equal shapes")
        effective_r = centres / (density * temperature)
        gas_method = "R=p/(rho*T) from simultaneous interpolated MCD fields"
    areoid = None
    if areoid_m is not None:
        areoid = _positive_array(areoid_m, "areoid radius", 2)
        if areoid.shape != orography.shape:
            raise ValueError("areoid radius and orography must have equal shapes")
        if np.any(areoid < 1.0e6) or np.any(areoid + orography <= 0.0):
            raise ValueError("areoid must be a planetary radius in metres, not an elevation offset")
    if geometry == "inverse_square" and areoid is None:
        raise ValueError("inverse_square geometry requires an explicit local areoid radius")
    interface_height = np.empty_like(interfaces)
    if geometry == "inverse_square":
        interface_height[0] = orography * areoid / (areoid + orography)
    else:
        interface_height[0] = orography
    height = np.empty_like(centres)
    scale_height = effective_r * temperature / gravity
    for level in range(centres.shape[0]):
        height[level] = interface_height[level] + scale_height[level] * np.log(interfaces[level] / centres[level])
        layer_log_pressure = np.log(interfaces[level] / interfaces[level + 1])
        interface_height[level + 1] = interface_height[level] + scale_height[level] * layer_log_pressure
    if geometry == "inverse_square":
        if np.any(interface_height >= areoid[None, :, :]):
            raise ValueError("hydrostatic column exceeds the spherical finite-potential domain")
        geometric = height * areoid[None, :, :] / (areoid[None, :, :] - height)
        geometric_interface = interface_height * areoid[None, :, :] / (areoid[None, :, :] - interface_height)
    else:
        geometric = height.copy()
        geometric_interface = interface_height.copy()
    reference_height = None if areoid is None else geometric + areoid[None, :, :] - reference_radius_m
    method = f"piecewise-isothermal hydrostatic; {geometry}; {gas_method}; terrain=source orography"
    return HeightResult(height, geometric, reference_height, interface_height,
        geometric_interface, effective_r, method)


def pressure_to_eta(pressure_pa: np.ndarray,
    surface_pressure_pa: np.ndarray, top_pressure_pa: float) -> np.ndarray:
    """Return (p-ptop)/(ps-ptop), a diagnostic, NOT real.exe's dry eta state."""
    pressure = _positive_array(pressure_pa, "pressure", 3)
    surface = _positive_array(surface_pressure_pa, "surface pressure", 2)
    top = float(top_pressure_pa)
    if pressure.shape[1:] != surface.shape:
        raise ValueError("pressure and surface pressure have incompatible horizontal shapes")
    if not np.isfinite(top) or top <= 0.0 or np.any(surface <= top):
        raise ValueError("top pressure must be positive and less than every surface pressure")
    if np.any(pressure > surface[None, :, :]) or np.any(pressure < top):
        raise ValueError("pressure lies outside the requested eta column")
    return (pressure - top) / (surface[None, :, :] - top)


def eta_interfaces_to_pressure(eta_interfaces: np.ndarray,
    surface_pressure_pa: np.ndarray, top_pressure_pa: float) -> np.ndarray:
    """Map explicit diagnostic eta interfaces 1..0 to pressure, without guessing WRF levels."""
    eta = _finite_array(eta_interfaces, "eta interfaces", 1)
    surface = _positive_array(surface_pressure_pa, "surface pressure", 2)
    top = float(top_pressure_pa)
    if eta.size < 2 or np.any(np.diff(eta) >= 0.0):
        raise ValueError("eta interfaces must be strictly decreasing from one to zero")
    if eta[0] != 1.0 or eta[-1] != 0.0 or np.any(eta < 0.0) or np.any(eta > 1.0):
        raise ValueError("explicit eta interfaces must start at one and finish at zero")
    if not np.isfinite(top) or top <= 0.0 or np.any(surface <= top):
        raise ValueError("top pressure must be positive and below every surface pressure")
    return top + eta[:, None, None] * (surface[None, :, :] - top)


def conservative_vertical_remap(mixing_ratio: np.ndarray,
    source_interfaces_pa: np.ndarray, target_interfaces_pa: np.ndarray) -> np.ndarray:
    """Pressure-overlap remap of nonnegative layer-mean mass mixing ratio.

    The result conserves integral(q*dp/g) over the shared target pressure range,
    under the declared piecewise-constant layer approximation and common g.
    Outside-domain extrapolation is forbidden; an explicitly truncated target
    preserves only mass within its pressure range, not the entire source column.
    """
    source = _descending_pressure(source_interfaces_pa, "source interfaces")
    target = _descending_pressure(target_interfaces_pa, "target interfaces")
    values = _finite_array(mixing_ratio, "mixing ratio", 3)
    if values.shape != (source.shape[0] - 1,) + source.shape[1:]:
        raise ValueError("mixing ratio and source interfaces describe different layers")
    if target.shape[1:] != source.shape[1:]:
        raise ValueError("source and target interfaces have incompatible horizontal shapes")
    if np.any(values < 0.0):
        raise ValueError("mass mixing ratio must be nonnegative")
    tolerance = np.maximum(source[0] * 1.0e-12, 1.0e-12)
    if np.any(target[0] > source[0] + tolerance) or np.any(target[-1] < source[-1] - tolerance):
        raise ValueError("target pressure column lies outside source; extrapolation is forbidden")
    output = np.empty((target.shape[0] - 1,) + target.shape[1:], dtype=np.float64)
    for level in range(output.shape[0]):
        lower = np.minimum(source[:-1], target[level][None, :, :])
        upper = np.maximum(source[1:], target[level + 1][None, :, :])
        overlap = np.maximum(lower - upper, 0.0)
        output[level] = np.sum(values * overlap, axis=0) / (target[level] - target[level + 1])
    return output


def interpolate_to_pressure(field: np.ndarray,
    pressure_center_pa: np.ndarray, target_pressure_pa: np.ndarray,
    surface_pressure_pa: np.ndarray | None = None, extrapolation: str = "missing") -> np.ndarray:
    """Log-pressure point interpolation; missing is default outside centre range.

    endpoint extrapolation must be explicitly requested, and remains forbidden
    below the supplied surface. This point interpolation must not replace the
    pressure-overlap remapper for mass-conserving aerosol transport.
    """
    centres = _descending_pressure(pressure_center_pa, "pressure centres")
    values = _finite_array(field, "field", 3)
    target = _positive_array(target_pressure_pa, "target pressures", 1)
    if values.shape != centres.shape:
        raise ValueError("field and pressure must have equal shapes")
    if extrapolation not in ("missing", "endpoint"):
        raise ValueError("extrapolation must be missing or endpoint")
    surface = None
    if surface_pressure_pa is not None:
        surface = _positive_array(surface_pressure_pa, "surface pressure", 2)
        if surface.shape != centres.shape[1:]:
            raise ValueError("surface pressure and pressure centres have incompatible horizontal shapes")
    output = np.empty((target.size,) + centres.shape[1:], dtype=np.float64)
    for latitude in range(centres.shape[1]):
        for longitude in range(centres.shape[2]):
            column = values[:, latitude, longitude][::-1]
            logarithm = np.log(centres[:, latitude, longitude][::-1])
            left = np.nan if extrapolation == "missing" else column[0]
            right = np.nan if extrapolation == "missing" else column[-1]
            interpolated = np.interp(np.log(target), logarithm, column, left=left, right=right)
            if surface is not None:
                interpolated[target > surface[latitude, longitude]] = np.nan
            output[:, latitude, longitude] = interpolated
    return output
