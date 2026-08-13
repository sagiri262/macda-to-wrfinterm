"""Mass-conserving MACDA dust initialization for two MarsWRF bins."""

from __future__ import annotations

import math

import numpy as np


def tracerFromDustMass(dustMixingRatio: np.ndarray, massFraction: float,
    tracerWeight: float) -> np.ndarray:
    """Return one tracer whose weighted value is its dust-bin mass."""

    return np.asarray(dustMixingRatio) * massFraction / tracerWeight


def combineDustTracers(tracer1: np.ndarray, tracer2: np.ndarray, tracerWeight1: float,
    tracerWeight2: float) -> np.ndarray:
    """
    重建 MarsWRF 的混合沙尘质量计算方式，0.3 * TRC01 + 0.7 * TRC02
    """

    return tracerWeight1 * np.asarray(tracer1) + tracerWeight2 * np.asarray(tracer2)


def validateDustPartition(dustMixingRatio: np.ndarray, tracer1: np.ndarray, tracer2: np.ndarray,
    tracerWeight1: float, tracerWeight2: float) -> None:
    """Reject a two-bin split that changes the total dust mixing ratio."""

    source = np.asarray(dustMixingRatio)
    combined = combineDustTracers(tracer1, tracer2, tracerWeight1, tracerWeight2)
    if not np.allclose(combined, source, rtol=1.0e-6, atol=1.0e-15):
        maximumError = float(np.max(np.abs(combined - source)))
        raise ValueError(f"Dust-bin partition does not conserve mass; "
                         f"maximum error={maximumError:.8g}")


def massExtinctionCoefficient(extinctionEfficiency: float, particleDensity: float,
    particleRadius: float) -> float:
    """Return spherical-particle mass extinction in square metres per kilogram."""

    return 0.75 * extinctionEfficiency / (particleDensity * particleRadius)


def opticalDepthScaling(macdaDensity: float, macdaEffectiveRadius: float, macdaExtinctionEfficiency: float,
    macdaGravity: float, marswrfRadiationDensity: float, marswrfRadiationGravity: float,
    dustDiameter1: float, dustDiameter2: float, extinctionEfficiency1: float,
    extinctionEfficiency2: float, massFraction1: float, massFraction2: float) -> float:
    """Return MarsWRF/MACDA opacity for an unchanged pressure-weighted dust column."""

    macdaCoefficient = massExtinctionCoefficient(macdaExtinctionEfficiency, macdaDensity,
        macdaEffectiveRadius)
    bin1Coefficient = massExtinctionCoefficient(extinctionEfficiency1, marswrfRadiationDensity,
        0.5 * dustDiameter1)
    bin2Coefficient = massExtinctionCoefficient(extinctionEfficiency2, marswrfRadiationDensity,
        0.5 * dustDiameter2)
    marswrfCoefficient = massFraction1 * bin1Coefficient + massFraction2 * bin2Coefficient
    return marswrfCoefficient * macdaGravity / (macdaCoefficient * marswrfRadiationGravity)


def sigmaLayerThickness(sigmaLevels: np.ndarray) -> np.ndarray:
    """Return sigma thicknesses for midpoint levels ordered surface to top."""

    levels = np.asarray(sigmaLevels, dtype=np.float64)
    if levels.ndim != 1 or levels.size < 2:
        raise ValueError("Sigma levels must be a one-dimensional array with at least two values")
    if np.any(levels <= 0.0) or np.any(levels > 1.0):
        raise ValueError("Sigma levels must be in the interval (0, 1]")
    if not np.all(np.diff(levels) < 0.0):
        raise ValueError("Sigma levels must be strictly descending from surface to top")

    boundaries = np.empty(levels.size + 1, dtype=np.float64)
    boundaries[0] = 1.0
    boundaries[-1] = 0.0
    boundaries[1:-1] = 0.5 * (levels[:-1] + levels[1:])
    return boundaries[:-1] - boundaries[1:]


def columnOpticalDepth(dustMixingRatio: np.ndarray, surfacePressure: np.ndarray, sigmaLevels: np.ndarray,
    massExtinction: float, gravity: float) -> np.ndarray:
    """Integrate dust optical depth on MACDA midpoint sigma levels."""

    dust = np.asarray(dustMixingRatio, dtype=np.float64)
    pressure = np.asarray(surfacePressure, dtype=np.float64)
    thickness = sigmaLayerThickness(sigmaLevels)
    if dust.shape[0] != thickness.size or dust.shape[1:] != pressure.shape:
        raise ValueError("Dust, surface pressure, and sigma dimensions are inconsistent")
    columnMass = pressure * np.sum(dust * thickness[:, None, None], axis=0) / gravity
    return massExtinction * columnMass


def validateDustParameters(massFraction1: float, massFraction2: float, tracerWeight1: float,
    tracerWeight2: float, positiveParameters: tuple[float, ...]) -> None:
    """Validate fractions, MarsWRF weights, and optical parameters."""

    fractions = (massFraction1, massFraction2)
    weights = (tracerWeight1, tracerWeight2)
    if not all(math.isfinite(value) for value in fractions + weights + positiveParameters):
        raise ValueError("Dust-bin fractions, weights, and optical parameters must be finite")
    if any(value < 0.0 for value in fractions):
        raise ValueError("Dust-bin mass fractions cannot be negative")
    if any(value <= 0.0 for value in weights + positiveParameters):
        raise ValueError("Dust-bin weights and optical parameters must be positive")
    if not math.isclose(sum(fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError("Dust-bin mass fractions must sum to 1.0")
    if not math.isclose(sum(weights), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError("MarsWRF dust-bin weights must sum to 1.0")
