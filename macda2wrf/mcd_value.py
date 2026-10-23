"""Signed MCD perturbations and optional coherent PCM EOF variability.

The default background remains the seasonally/diurnally interpolated mean.
RMS is an unsigned scale, not a daily anomaly. Optional EOF samples represent
PCM weather variability, not observations of the actual MY36 weather.
"""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path

import numpy as np


MCD_ORBIT_SOLS = 668.5921
EOF_VARIANCE_FACTOR = 1.205
EOF_TIME_WINDOW_SOLS = 30.0
REFERENCE_PRESSURE_PA = 610.0
ANOMALY_CONVENTIONS = ("signed_anomaly", "estimate_minus_truth_bias")
# 温度、UVW风速三分量、气压、密度、地表温度、沙尘光学厚度
RMS_UNITS = {"temp": "K",
             "u": "m s-1",
             "v": "m s-1",
             "w": "m s-1",
             "ps": "Pa",
             "rho": "kg m-3",
             "tsurf": "K",
             "tau_pref_gcm": "1"}


def apply_anomaly(mean: np.ndarray, anomaly: np.ndarray,
    convention: str) -> np.ndarray:
    """Combine fields only with an explicit, documented sign definition.

    ``signed_anomaly`` means anomaly = target - mean, hence target = mean +
    anomaly. ``estimate_minus_truth_bias`` means bias = estimate - target,
    hence target = estimate - bias. Neither option makes the target a truth
    measurement. Inputs must have equal shapes; RMS is never accepted as a
    convention or silently assigned a sign.
    """
    if convention not in ANOMALY_CONVENTIONS:
        raise ValueError("convention must be signed_anomaly or estimate_minus_truth_bias; RMS is unsigned")
    background = np.asarray(mean, dtype=np.float64)
    perturbation = np.asarray(anomaly, dtype=np.float64)
    if background.shape != perturbation.shape:
        raise ValueError("mean and anomaly must have identical shapes")
    if not np.all(np.isfinite(background)) or not np.all(np.isfinite(perturbation)):
        raise ValueError("mean and anomaly must contain finite values")
    if convention == "signed_anomaly":
        return background + perturbation
    return background - perturbation


def _coordinate(values: np.ndarray, name: str) -> np.ndarray:
    # 坐标的数据类型
    coordinate = np.asarray(values, dtype=np.float64)
    # 弹出错误
    if coordinate.ndim != 1 or coordinate.size == 0 or not np.all(np.isfinite(coordinate)):
        raise ValueError(f"{name} must be a finite nonempty one-dimensional coordinate")
    return coordinate


def read_rms_fields(rms_path: str | Path, variables: tuple = ("temp", "u", "v", "ps", "tau_pref_gcm")) -> dict:
    """Read nonnegative native-grid weather scales, never assign anomaly signs.

    RMS = sqrt(mean((one-sol mean - ten-sol mean)**2)) in the MCD reference
    definition. These scales do not identify the actual day-specific anomaly,
    independent truth error, or off-diagonal elements of B. Empty file units
    are annotated using the matching physical quantity's documented units.
    """
    from scipy.io import netcdf_file

    path = Path(rms_path).expanduser()
    if not path.is_file():
        raise FileNotFoundError(f"RMS file not found: {path}")
    requested = tuple(variables)
    if not requested or any(name not in RMS_UNITS for name in requested):
        raise ValueError("request documented rms variables; no rmsdustq is supplied by this interface")
    fields = {}
    with netcdf_file(str(path), mmap=True) as dataset:
        for name in requested:
            source_name = "rms" + name
            if source_name not in dataset.variables:
                raise ValueError(f"RMS file lacks {source_name}")
            values = np.array(dataset.variables[source_name][:], dtype=np.float64, copy=True)
            if values.ndim not in (2, 3) or not np.all(np.isfinite(values)) or np.any(values < 0.0):
                raise ValueError(f"{source_name} must be a finite nonnegative horizontal field or profile")
            fields[name] = values
    fields["metadata"] = {"source_file": str(path.resolve()), "kind": "unsigned_pcm_weather_variability_scale",
        "definition": "sqrt(mean((one_sol_mean-ten_sol_mean)**2))", "units": dict(RMS_UNITS),
        "limitations": "not signed daily anomalies, not a validated background-error covariance"}
    return fields


def _linear_indices(source: np.ndarray, target: np.ndarray) -> tuple:
    """Return bracketing indices, clipping outside the source latitude/height."""
    if source.size == 1:
        zero = np.zeros(target.shape, dtype=np.int64)
        return zero, zero, np.zeros(target.shape)
    upper = np.clip(np.searchsorted(source, target, side="right"), 1, source.size - 1)
    lower = upper - 1
    fraction = (target - source[lower]) / (source[upper] - source[lower])
    return lower, upper, np.clip(fraction, 0.0, 1.0)


class McdVariability:
    """Read and sample additive ``1.205 * (pc-pcsmth) * EOF * norm``.

    One seeded phase shift within +/-15 sols is fixed for the whole instance,
    shared by all variables and grid points. No independent grid-cell noise
    is introduced. The EOF coordinate called altitude is used only to map
    modes between native PCM levels; it is not a geometric/geopotential height.
    Results have shapes (level, latitude, longitude), or (latitude, longitude)
    for ps. No dust perturbation or dust-weather covariance is provided.
    """

    def __init__(self, eof_path: str | Path, seed: int = 0, amplitude: float = 1.0, modes: int = 200, orbit_sols: float = MCD_ORBIT_SOLS) -> None:
        self.eof_path = Path(eof_path).expanduser()
        if not self.eof_path.is_file():
            raise FileNotFoundError(f"EOF file not found: {self.eof_path}")
        if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
            raise ValueError("seed must be a nonnegative integer")
        if isinstance(modes, bool) or not isinstance(modes, (int, np.integer)) or modes < 1:
            raise ValueError("modes must be a positive integer")
        if not np.isfinite(amplitude) or amplitude < 0.0:
            raise ValueError("amplitude must be finite and nonnegative")
        if not np.isfinite(orbit_sols) or orbit_sols <= 1.0:
            raise ValueError("orbit_sols must be finite and greater than one")
        self.seed = int(seed)
        self.amplitude = float(amplitude)
        self.requested_modes = int(modes)
        self.orbit_sols = float(orbit_sols)
        generator = np.random.default_rng(self.seed)
        self.phase_shift_sols = float((generator.random() - 0.5) * EOF_TIME_WINDOW_SOLS)
        self._loaded = False
        self._daily_cache = OrderedDict()

    def _load(self) -> None:
        if self._loaded:
            return
        from scipy.io import netcdf_file

        names = ("latitude", "longitude", "altitude", "time", "pc", "pcsmth")
        names += ("temp", "u", "v", "ps", "normt", "normu", "normv", "normp")
        with netcdf_file(str(self.eof_path), mmap=True) as dataset:
            missing = [name for name in names if name not in dataset.variables]
            if missing:
                raise ValueError(f"EOF file lacks required variables: {', '.join(missing)}")
            data = {name: np.array(dataset.variables[name][:], copy=True) for name in names}
        self.latitude = _coordinate(data["latitude"], "EOF latitude")
        self.longitude = _coordinate(data["longitude"], "EOF longitude")
        self.altitude = _coordinate(data["altitude"], "EOF pseudo-altitude")
        self.times = _coordinate(data["time"], "EOF time")
        if self.latitude.size < 2 or self.longitude.size < 2 or self.altitude.size < 2:
            raise ValueError("EOF spatial coordinates need at least two points per axis")
        if np.any(np.diff(self.altitude) <= 0.0):
            raise ValueError("EOF pseudo-altitude must increase surface to top")
        if np.any(np.diff(np.sort(self.latitude)) <= 0.0):
            raise ValueError("EOF latitude must not contain duplicate values")
        if np.any(np.abs(self.latitude) > 90.0):
            raise ValueError("EOF latitude lies outside [-90, 90]")
        if not np.allclose(self.times, np.arange(1, self.times.size + 1)):
            raise ValueError("EOF time must contain consecutive model days starting at one")
        if self.times.size != int(np.ceil(self.orbit_sols)):
            raise ValueError("EOF last day count does not match the stated model orbital year")
        mode_count = data["pc"].shape[0]
        self.mode_count = min(self.requested_modes, mode_count)
        expected_pc = (mode_count, self.times.size, self.latitude.size)
        if data["pc"].shape != expected_pc or data["pcsmth"].shape != expected_pc:
            raise ValueError("EOF principal components must have dimensions (mode, day, latitude)")
        self._latitude_order = np.argsort(self.latitude)
        wrapped_longitude = self.longitude % 360.0
        self._longitude_order = np.argsort(wrapped_longitude)
        self._longitude_sorted = wrapped_longitude[self._longitude_order]
        if np.any(np.diff(self._longitude_sorted) <= 0.0):
            raise ValueError("EOF longitude must not duplicate the periodic seam")
        self._latitude_sorted = self.latitude[self._latitude_order]
        selected = slice(0, self.mode_count)
        self._principal = EOF_VARIANCE_FACTOR * (data["pc"][selected] - data["pcsmth"][selected])
        self._principal = np.asarray(self._principal, dtype=np.float32)
        self._modes = {}
        self._normalization = {}
        norm_names = {"temp": "normt", "u": "normu", "v": "normv", "ps": "normp"}
        for name, norm_name in norm_names.items():
            expected_shape = (mode_count, self.latitude.size, self.longitude.size)
            if name != "ps":
                expected_shape = (mode_count, self.altitude.size, self.latitude.size, self.longitude.size)
            if data[name].shape != expected_shape or data[norm_name].shape != (self.latitude.size,):
                raise ValueError(f"EOF dimensions are inconsistent for {name}")
            self._modes[name] = np.asarray(data[name][selected], dtype=np.float32)
            self._normalization[name] = np.asarray(data[norm_name], dtype=np.float32)
            if not np.all(np.isfinite(self._modes[name])):
                raise ValueError(f"EOF modes contain nonfinite values for {name}")
            if not np.all(np.isfinite(self._normalization[name])):
                raise ValueError(f"EOF normalization contains nonfinite values for {name}")
        if not np.all(np.isfinite(self._principal)):
            raise ValueError("EOF principal components contain nonfinite values")
        self._loaded = True

    def _time_indices(self, elapsed_sols: float) -> tuple:
        phase = (elapsed_sols + self.phase_shift_sols) % self.orbit_sols
        last_integer = self.times.size - 1
        if phase < 1.0:
            return self.times.size - 1, 0, phase, phase
        if phase >= last_integer:
            weight = (phase - last_integer) / (self.orbit_sols - last_integer)
            return last_integer - 1, last_integer, weight, phase
        lower_day = int(np.floor(phase))
        return lower_day - 1, lower_day, phase - lower_day, phase

    def _native_day(self, index: int) -> dict:
        if index in self._daily_cache:
            self._daily_cache.move_to_end(index)
            return self._daily_cache[index]
        principal = self._principal[:, index, :]
        fields = {}
        for name, eigenvectors in self._modes.items():
            if name == "ps":
                field = np.einsum("ml,mly->ly", principal, eigenvectors, optimize=True)
                field = field * self._normalization[name][:, None]
            else:
                field = np.einsum("ml,mzly->zly", principal, eigenvectors, optimize=True)
                field = field * self._normalization[name][None, :, None]
            fields[name] = field[..., self._latitude_order, :][..., self._longitude_order]
        self._daily_cache[index] = fields
        while len(self._daily_cache) > 3:
            self._daily_cache.popitem(last=False)
        return fields

    def _horizontal(self, field: np.ndarray, latitude: np.ndarray,
        longitude: np.ndarray) -> np.ndarray:
        lower_lat, upper_lat, weight_lat = _linear_indices(self._latitude_sorted, latitude)
        latitude_field = (1.0 - weight_lat[:, None]) * field[..., lower_lat, :]
        latitude_field += weight_lat[:, None] * field[..., upper_lat, :]
        longitude_extended = np.r_[self._longitude_sorted, self._longitude_sorted[0] + 360.0]
        periodic_field = np.concatenate((latitude_field, latitude_field[..., :1]), axis=-1)
        target_longitude = (longitude - self._longitude_sorted[0]) % 360.0 + self._longitude_sorted[0]
        lower_lon, upper_lon, weight_lon = _linear_indices(longitude_extended, target_longitude)
        result = (1.0 - weight_lon) * periodic_field[..., lower_lon]
        return result + weight_lon * periodic_field[..., upper_lon]

    def sample(self, elapsed_sols: float, latitude: np.ndarray,
        longitude: np.ndarray, aps: np.ndarray, bps: np.ndarray, source_ps: np.ndarray | None = None) -> dict:
        """Return additive PCM variability and explicit provenance metadata.

        ``aps`` and ``bps`` identify target hybrid levels through MCD's
        -10*log(aps/610+bps) pseudo-coordinate. Surface pressure is only checked
        for compatible shape/positivity: it does not turn this coordinate into
        actual height. Polar values outside the EOF latitude range use its
        nearest latitude, matching official spatial boundary handling.
        """
        if not np.isfinite(elapsed_sols):
            raise ValueError("elapsed_sols must be finite")
        self._load()
        target_latitude = _coordinate(latitude, "target latitude")
        target_longitude = _coordinate(longitude, "target longitude")
        pressure_term = _coordinate(aps, "aps")
        sigma_term = _coordinate(bps, "bps")
        if pressure_term.shape != sigma_term.shape:
            raise ValueError("aps and bps must have equal lengths")
        if np.any(np.abs(target_latitude) > 90.0):
            raise ValueError("target latitude lies outside [-90, 90]")
        pseudo_pressure = pressure_term / REFERENCE_PRESSURE_PA + sigma_term
        if np.any(pseudo_pressure <= 0.0):
            raise ValueError("aps/610+bps must be strictly positive")
        target_altitude = -10.0 * np.log(pseudo_pressure)
        if np.any(np.diff(target_altitude) <= 0.0):
            raise ValueError("target hybrid levels must increase in pseudo-altitude")
        if source_ps is not None:
            surface = np.asarray(source_ps, dtype=np.float64)
            expected_shape = (target_latitude.size, target_longitude.size)
            if surface.shape != expected_shape or not np.all(np.isfinite(surface)) or np.any(surface <= 0.0):
                raise ValueError("source_ps must be finite positive pressure on the requested horizontal grid")
        lower, upper, time_weight, phase = self._time_indices(float(elapsed_sols))
        lower_fields = self._native_day(lower)
        upper_fields = self._native_day(upper)
        lower_z, upper_z, weight_z = _linear_indices(self.altitude, target_altitude)
        result = {}
        for name in ("temp", "u", "v", "ps"):
            field = (1.0 - time_weight) * lower_fields[name] + time_weight * upper_fields[name]
            field = self._horizontal(field, target_latitude, target_longitude)
            if name != "ps":
                interpolated = (1.0 - weight_z[:, None, None]) * field[lower_z]
                field = interpolated + weight_z[:, None, None] * field[upper_z]
            result[name] = self.amplitude * field
        metadata = {"kind": "synthetic_pcm_eof_variability", "source_file": str(self.eof_path.resolve()),
            "formula": "amplitude*1.205*sum((pc-pcsmth)*EOF)*latitude_norm", "seed": self.seed,
            "phase_shift_sols": self.phase_shift_sols, "sample_phase_sols": phase,
            "amplitude": self.amplitude, "mode_count": self.mode_count, "requested_modes": self.requested_modes,
            "sign_convention": "signed_anomaly: target=mean+anomaly", "dust_cross_covariance": "unavailable",
            "height_mapping": "pseudo-altitude for EOF mode interpolation only; not GHT",
            "time_boundary": f"last EOF sample at orbit end; final fractional sol normalized to {self.orbit_sols}",
            "pole_boundary": "nearest EOF latitude; native EOF excludes exact poles",
            "truth_claim": "PCM variability sample, not recovered observed day-specific weather"}
        result["metadata"] = metadata
        return result
