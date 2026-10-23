"""Terminal-speed diagnostics, not a sedimentation time integration.

MarsWRF ``stokes_dust`` consumes particle DIAMETER and calculates a positive
downward speed. Its two bins are discrete sizes with independently transported,
unweighted mass tracers. The reference lognormal calculation below is an
explicit distribution assumption: an effective radius alone does not determine
the mass-weighted terminal speed.

Source: WRFV3/phys/module_mp_sedim_dust.F:186-242 and :147-162;
WRFV3/phys/module_microphysics_driver.F:1499-1509;
MCD_6.2/mcd/MCD.F90:6470-6525. No radiation weights or floors enter here.
"""

from __future__ import annotations

import numpy as np


CO2_MOLECULAR_RADIUS_M = 2.2e-10
WRF_MEAN_FREE_PATH_COEFFICIENT = 0.707 * 8.31 / (4.0 * 3.1416 * CO2_MOLECULAR_RADIUS_M**2 * 6.023e23)
WRF_SLIP_RADIUS_COEFFICIENT = 2.6667 / 2.0
MCD_REFERENCE_GRAVITY_M_S2 = 3.7257964


def _positive_input(value: np.ndarray, name: str) -> np.ndarray:
    """Retain NaNs as missing; reject infinity and nonpositive physical data."""

    result = np.asarray(value, dtype=np.float64)
    if np.any(np.isinf(result)) or np.any(result <= 0.0):
        raise ValueError(f"{name} must be positive or NaN")
    return result


def _positive_constant(value: float, name: str) -> float:
    result = np.asarray(value, dtype=np.float64)
    if result.ndim != 0 or not np.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be a finite positive scalar")
    return float(result)


def _coefficients(particle_density: float,
    viscosity: float, gravity_m_s2: float) -> float:
    density = _positive_constant(particle_density, "particle_density")
    viscosity = _positive_constant(viscosity, "viscosity")
    gravity = _positive_constant(gravity_m_s2, "gravity_m_s2")
    return (2.0 / 9.0) * density * gravity / viscosity


def wrf_stokes_velocity(pressure_pa: np.ndarray,
    temperature_k: np.ndarray, radius_m: np.ndarray, particle_density: float = 2500.0,
    viscosity: float = 1.0e-5, gravity_m_s2: float = 3.727) -> np.ndarray:
    """Reproduce ``stokes_dust`` algebra in float64, returning downward m/s.

    ``p``, ``T`` and radius broadcast; they must be positive or NaN. The radius
    is converted to particle diameter before evaluating the source law. The
    molecular-viscosity input is fixed, not temperature-dependent. All scalar
    constants must be finite and positive. The source caches its coefficient
    on first invocation; this stateless diagnostic evaluates each explicitly
    supplied configuration independently instead of reproducing that cache.

    The model calls this at LOWER-interface pressure ``p8w(k)`` with layer
    temperature ``theta(k)*pi(k)``; model-like evaluation must pair those
    inputs deliberately. A layer-center-pressure comparison is a separately
    labelled sensitivity, not the exact model calling convention.
    """

    pressure, temperature, radius = np.broadcast_arrays(_positive_input(pressure_pa, "pressure_pa"),
        _positive_input(temperature_k, "temperature_k"), _positive_input(radius_m, "radius_m"))
    diameter = 2.0 * radius
    coefficient = _coefficients(particle_density, viscosity, gravity_m_s2)
    mean_free_path = WRF_MEAN_FREE_PATH_COEFFICIENT * temperature / pressure
    return 0.25 * coefficient * diameter**2 * (1.0 + 2.6667 * mean_free_path / diameter)


def wrf_lognormal_mass_velocity(pressure_pa: np.ndarray,
    temperature_k: np.ndarray, reff_m: np.ndarray, effective_variance: float = 0.5,
    particle_density: float = 2500.0, viscosity: float = 1.0e-5,
    gravity_m_s2: float = 3.727) -> np.ndarray:
    """Exact lognormal mass mean for the WRF spherical-particle speed law.

    This is a chosen distribution closure, NOT an observed reference speed.
    With number moments ``M_k``, ``reff=M3/M2`` and
    ``v_eff=M4*M2/M3**2-1``, lognormal ``exp(sigma_log**2)=1+v_eff``.
    Consequently ``M5/M3=reff**2*(1+v_eff)**3`` and
    ``M4/M3=reff*(1+v_eff)``. They weight the quadratic Stokes and linear slip
    terms separately. Setting ``v_eff=0`` gives a single-size comparison at
    ``reff``; using ``v(reff)`` for an unspecified distribution is not valid.
    """

    variance = np.asarray(effective_variance, dtype=np.float64)
    if variance.ndim != 0 or not np.isfinite(variance) or variance < 0.0:
        raise ValueError("effective_variance must be a finite nonnegative scalar")
    pressure, temperature, radius = np.broadcast_arrays(_positive_input(pressure_pa, "pressure_pa"),
        _positive_input(temperature_k, "temperature_k"), _positive_input(reff_m, "reff_m"))
    coefficient = _coefficients(particle_density, viscosity, gravity_m_s2)
    width = 1.0 + float(variance)
    mean_free_path = WRF_MEAN_FREE_PATH_COEFFICIENT * temperature / pressure
    return coefficient * (radius**2 * width**3 + WRF_SLIP_RADIUS_COEFFICIENT * mean_free_path * radius * width)


def mcd_deposition_velocity(surface_pressure_pa: np.ndarray,
    temperature_k: np.ndarray, reff_m: np.ndarray, gravity_m_s2: float = 3.727) -> np.ndarray:
    """MCD ``dust_deposition`` speed recipe only, returning downward m/s.

    The official diagnostic is a NEAR-SURFACE deposition approximation. It
    fixes density=2500 kg/m3, viscosity=1e-5 Pa s, ``v_eff=0.5`` and
    ``r_sed=(1+v_eff)**2*reff``; its molecular coefficient is multiplied by
    0.5 for a disk-shape correction and its slip coefficient is 1.333.
    This recipe is NOT the exact lognormal mass-mean speed above, and extending
    it aloft is only a recipe sensitivity, not a validated MCD profile output.
    MCD deposition itself is ``mmr*rho_air*speed``; it is not computed here.
    Default gravity matches the WRF diagnostic for controlled comparison;
    pass ``MCD_REFERENCE_GRAVITY_M_S2`` to reproduce the native MCD constant.
    """

    pressure, temperature, radius = np.broadcast_arrays(_positive_input(surface_pressure_pa, "surface_pressure_pa"),
        _positive_input(temperature_k, "temperature_k"), _positive_input(reff_m, "reff_m"))
    sedimentation_radius = (1.0 + 0.5)**2 * radius
    gravity = _positive_constant(gravity_m_s2, "gravity_m_s2")
    coefficient = (2.0 / 9.0) * 2500.0 * gravity / 1.0e-5
    corrected_path = 0.5 * WRF_MEAN_FREE_PATH_COEFFICIENT * temperature / pressure
    return coefficient * sedimentation_radius**2 * (1.0 + 1.333 * corrected_path / sedimentation_radius)
