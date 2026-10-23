"""Independently validate MCD background NetCDF files without writing reports.

The validator samples full three-dimensional records, but reads every time
coordinate. It checks the discrete column represented by the inferred pressure
interfaces; passing does not establish native PCM finite-volume bounds, actual
Mars weather, WRF dry-eta initial conditions, or a completed reanalysis.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from typing import Any

import numpy as np


REQUIRED_DIMENSIONS = ("Time", "level", "interface", "latitude", "longitude", "DateStrLen")
VOLUME_VARIABLES = ("pressure", "eta_diagnostic", "ght", "temp", "u", "v", "rho", "dustq", "reffdust", "trc01", "trc02")
COLUMN_VARIABLES = ("ps", "dust_column_mass", "tau_visible610_mcd", "tau_visible610_spherical", "tau_wrf_native700")
TIME_VARIABLES = ("Time", "solar_longitude", "mean_solar_time", "true_solar_time")
UNIT_ALIASES = {"pressure": ("Pa",), "pressure_interface": ("Pa",), "ps": ("Pa",), "ght": ("m",),
    "temp": ("K",), "u": ("m/s", "m s-1", "m s^-1"), "v": ("m/s", "m s-1", "m s^-1"),
    "rho": ("kg/m3", "kg m-3", "kg m^-3", "kg/m^3"),
    "dustq": ("kg/kg", "kg kg-1", "kg kg^-1", "1"), "reffdust": ("m",),
    "trc01": ("kg/kg", "kg kg-1", "kg kg^-1", "1"),
    "trc02": ("kg/kg", "kg kg-1", "kg kg^-1", "1"),
    "dust_column_mass": ("kg/m2", "kg m-2", "kg m^-2", "kg/m^2"),
    "orography": ("m",), "areoid_radius": ("m",), "eta_diagnostic": ("1",),
    "tau_visible610_mcd": ("1",), "tau_visible610_spherical": ("1",), "tau_wrf_native700": ("1",)}
DATE_PATTERN = re.compile(r"^(\d{4})-(\d{5})_(\d{2}):(\d{2}):(\d{2})$")
GRAVITY_M_S2 = 3.727
SAMPLES_PER_SOL = 12
MY36_SOL_COUNT = 669
SI_SECONDS_PER_SOL = 88775.245


class VariableView:
    """Expose either supported NetCDF schema through one read-only contract."""

    def __init__(self, source: Any, name: str,
        dimension_names: dict[str, str], si_time: bool = False) -> None:
        self.source = source
        self.name = name
        self.si_time = si_time
        self.dimensions = tuple(dimension_names.get(value, value) for value in source.dimensions)
        if name == "Times":
            self.dimensions = ("Time", "DateStrLen")

    def __getattr__(self, name: str) -> Any:
        return getattr(self.source, name)

    def __getitem__(self, index: Any) -> Any:
        values = self.source[index]
        if self.name == "Time" and self.si_time:
            return np.ma.asarray(values) / SI_SECONDS_PER_SOL
        if self.name == "Times" and np.ndim(values) == 1:
            return np.asarray([list(str(value)) for value in values], dtype="S1")
        return values


class DatasetView:
    """Support canonical legacy-style names and the pipeline's lowercase names."""

    def __init__(self, source: Any) -> None:
        self.source = source
        dimension_names = {"time": "Time", "lat": "latitude", "lon": "longitude"}
        aliases = {"time": "Time", "lat": "latitude", "lon": "longitude", "temperature": "temp",
            "q_bin1": "trc01", "q_bin2": "trc02", "column_dust_mass": "dust_column_mass",
            "ls": "solar_longitude", "model_hour": "mean_solar_time", "wrf_date": "Times",
            "ltst0": "true_solar_time"}
        self.dimensions = {dimension_names.get(name, name): dimension for name, dimension in source.dimensions.items()}
        self.variables = {aliases.get(name, name): VariableView(variable, aliases.get(name, name),
            dimension_names, si_time=name == "time") for name, variable in source.variables.items()}
        if "wrf_date" in source.variables:
            self.dimensions.setdefault("DateStrLen", range(19))

    def __getattr__(self, name: str) -> Any:
        return getattr(self.source, name)


def _add_error(result: dict[str, Any], message: str) -> None:
    result["errors"].append(message)


def _array(variable: Any, index: Any = slice(None)) -> np.ndarray:
    values = np.ma.asarray(variable[index])
    return np.asarray(np.ma.filled(values, np.nan), dtype=np.float64)


def _finite(result: dict[str, Any], name: str,
    values: np.ndarray, positive: bool = False) -> bool:
    if np.any(~np.isfinite(values)):
        _add_error(result, f"{name} contains missing or nonfinite values")
        return False
    if positive and np.any(values <= 0.0):
        _add_error(result, f"{name} must be strictly positive")
        return False
    return True


def _close(result: dict[str, Any], name: str,
    actual: np.ndarray, expected: np.ndarray, rtol: float = 3.0e-5, atol: float = 1.0e-10) -> None:
    if actual.shape != expected.shape or not np.allclose(actual, expected, rtol=rtol, atol=atol):
        _add_error(result, f"{name} does not close within rtol={rtol:g}, atol={atol:g}")


def _check_metadata(dataset: Any, result: dict[str, Any]) -> bool:
    missing = [name for name in REQUIRED_DIMENSIONS if name not in dataset.dimensions]
    if missing:
        _add_error(result, f"missing dimensions: {', '.join(missing)}")
        return False
    lengths = {name: len(dataset.dimensions[name]) for name in REQUIRED_DIMENSIONS}
    result["dimensions"] = lengths
    if lengths["interface"] != lengths["level"] + 1:
        _add_error(result, "interface dimension must equal level + 1")
    if lengths["Time"] < 1 or lengths["level"] < 2:
        _add_error(result, "at least one time and two vertical layers are required")
        return False
    names = TIME_VARIABLES[:-1] + VOLUME_VARIABLES + COLUMN_VARIABLES
    names += ("Times", "pressure_interface", "latitude", "longitude", "orography", "areoid_radius")
    absent = [name for name in names if name not in dataset.variables]
    if absent:
        _add_error(result, f"missing variables: {', '.join(absent)}")
        return False
    volume_dims = ("Time", "level", "latitude", "longitude")
    column_dims = ("Time", "latitude", "longitude")
    expected = {name: volume_dims for name in VOLUME_VARIABLES}
    expected.update({name: column_dims for name in COLUMN_VARIABLES})
    expected.update({name: ("Time",) for name in TIME_VARIABLES if name in dataset.variables})
    expected["pressure_interface"] = ("Time", "interface", "latitude", "longitude")
    expected["Times"] = ("Time", "DateStrLen")
    expected["orography"] = ("latitude", "longitude")
    expected["areoid_radius"] = ("latitude", "longitude")
    for name, dimensions in expected.items():
        if dataset.variables[name].dimensions != dimensions:
            _add_error(result, f"{name} dimensions must be {dimensions}")
    for name, aliases in UNIT_ALIASES.items():
        units = str(getattr(dataset.variables[name], "units", "")).strip()
        if units not in aliases:
            _add_error(result, f"{name} units {units!r} are not one of {aliases}")
    product = str(getattr(dataset, "product_type", ""))
    if not product:
        result["warnings"].append("product_type is absent; background/analysis provenance is unspecified")
    method = str(getattr(dataset, "value_method", "")).lower().replace(" ", "_")
    if any(value in method for value in ("mean_minus_rms", "mean-rms", "mean_subtract_rms")):
        _add_error(result, "mean minus RMS cannot reconstruct a dated true state")
    eta_text = " ".join(str(getattr(dataset.variables["eta_diagnostic"], key, ""))
        for key in ("long_name", "description", "comment", "diagnostic_only"))
    if hasattr(dataset.variables["eta_diagnostic"], "diagnostic_only"):
        eta_text += " diagnostic"
    if "diagnostic" not in eta_text.lower():
        result["warnings"].append("eta_diagnostic should explicitly state diagnostic-only status")
    gravity = float(getattr(dataset, "gravity_m_s2", GRAVITY_M_S2))
    if not np.isfinite(gravity) or gravity <= 0.0:
        _add_error(result, "gravity_m_s2 must be positive and finite")
    if not np.isclose(gravity, GRAVITY_M_S2, rtol=0.0, atol=1.0e-8):
        result["warnings"].append(f"declared gravity {gravity:g} differs from target MarsWRF 3.727")
    reference = float(getattr(dataset, "reference_pressure_pa", 610.0))
    if reference != 610.0:
        _add_error(result, "reference_pressure_pa must be 610 Pa; historical tau700 is a separate diagnostic")
    return not result["errors"]


def _check_time(dataset: Any, result: dict[str, Any]) -> np.ndarray:
    time_variable = dataset.variables["Time"]
    if time_variable.si_time and str(getattr(time_variable.source, "units", "")) not in ("s", "seconds"):
        _add_error(result, "lowercase time coordinate must declare SI seconds")
    times = _array(dataset.variables["Time"])
    if not _finite(result, "Time", times):
        return times
    if np.any(times < 0.0) or np.any(times >= MY36_SOL_COUNT):
        _add_error(result, "Time lies outside the MY36 model-label year [0, 669) sols")
    if times.size > 1 and not np.allclose(np.diff(times), 1.0 / SAMPLES_PER_SOL, rtol=0.0, atol=1.0e-7):
        _add_error(result, "Time must increase by exactly two labelled Martian hours (1/12 sol)")
    chars = np.asarray(dataset.variables["Times"][:])
    if chars.dtype.kind == "S":
        labels = [b"".join(row).decode("ascii").rstrip("\0 ") for row in chars]
    else:
        labels = ["".join(row).rstrip("\0 ") for row in chars.astype(str)]
    for position, (elapsed, label) in enumerate(zip(times, labels)):
        match = DATE_PATTERN.fullmatch(label)
        if match is None:
            _add_error(result, f"Times[{position}] is not YYYY-DDDDD_HH:MM:SS")
            break
        year, sol, hour, minute, second = map(int, match.groups())
        label_elapsed = (sol - 1) + (hour * 3600 + minute * 60 + second) / 86400.0
        if year != 36 or not 1 <= sol <= MY36_SOL_COUNT or hour >= 24 or minute >= 60 or second >= 60:
            _add_error(result, f"Times[{position}] is outside the MY36 labelled calendar")
            break
        if not np.isclose(elapsed, label_elapsed, rtol=0.0, atol=1.0e-7):
            _add_error(result, f"Times[{position}] and elapsed Time refer to different labelled instants")
            break
    for name in TIME_VARIABLES[1:]:
        if name not in dataset.variables:
            result["warnings"].append(f"{name} is absent; its range could not be checked")
            continue
        values = _array(dataset.variables[name])
        upper = 360.0 if name == "solar_longitude" else 24.0
        if _finite(result, name, values) and (np.any(values < 0.0) or np.any(values >= upper)):
            _add_error(result, f"{name} must be within [0, {upper:g})")
    return times


def _check_grid(dataset: Any, result: dict[str, Any]) -> None:
    for name, expected_spacing in (("latitude", 3.75), ("longitude", 5.625)):
        values = _array(dataset.variables[name])
        if not _finite(result, name, values):
            continue
        differences = np.diff(values)
        if values.size < 2 or not (np.all(differences > 0.0) or np.all(differences < 0.0)):
            _add_error(result, f"{name} must be strictly monotonic")
            continue
        if not np.allclose(np.abs(differences), expected_spacing):
            result["warnings"].append(f"{name} differs from native MCD {expected_spacing:g}-degree spacing")
    orography = _array(dataset.variables["orography"])
    radius = _array(dataset.variables["areoid_radius"])
    _finite(result, "orography", orography)
    if _finite(result, "areoid_radius", radius) and np.any(radius < 1.0e6):
        _add_error(result, "areoid_radius must be a planetary radius in metres, not an elevation offset")


def _check_record(dataset: Any, result: dict[str, Any], index: int) -> None:
    fields = {name: _array(dataset.variables[name], index) for name in VOLUME_VARIABLES + COLUMN_VARIABLES}
    fields["pressure_interface"] = _array(dataset.variables["pressure_interface"], index)
    for name, values in fields.items():
        positive = name in ("pressure", "pressure_interface", "ps", "temp", "rho")
        _finite(result, f"{name}[{index}]", values, positive)
    nonnegative = ("dustq", "trc01", "trc02", "dust_column_mass", "tau_visible610_mcd",
        "tau_visible610_spherical", "tau_wrf_native700")
    for name in nonnegative:
        if np.any(fields[name] < 0.0):
            _add_error(result, f"{name}[{index}] must be nonnegative")
    invalid_radius = (fields["reffdust"] < 0.0) | ((fields["dustq"] > 0.0) & (fields["reffdust"] <= 0.0))
    if np.any(invalid_radius):
        _add_error(result, f"reffdust[{index}] must be positive wherever dust mass is positive")
    pressure = fields["pressure"]
    interfaces = fields["pressure_interface"]
    if interfaces.shape[0] != pressure.shape[0] + 1:
        _add_error(result, "pressure centre/interface layer counts differ")
        return
    if np.any(np.diff(pressure, axis=0) >= 0.0) or np.any(np.diff(interfaces, axis=0) >= 0.0):
        _add_error(result, f"pressure at record {index} is not strictly decreasing surface to top")
    if np.any(pressure >= interfaces[:-1]) or np.any(pressure <= interfaces[1:]):
        _add_error(result, f"pressure centres at record {index} must lie strictly inside their interfaces")
    _close(result, f"surface interface[{index}]", interfaces[0], fields["ps"], atol=1.0e-4)
    expected_eta = (pressure - interfaces[-1]) / (fields["ps"] - interfaces[-1])
    _close(result, f"diagnostic eta[{index}]", fields["eta_diagnostic"], expected_eta, atol=1.0e-7)
    if np.any(fields["eta_diagnostic"] < 0.0) or np.any(fields["eta_diagnostic"] > 1.0):
        _add_error(result, f"eta_diagnostic[{index}] must lie in [0, 1]")
    height = fields["ght"]
    if np.any(np.diff(height, axis=0) <= 0.0):
        _add_error(result, f"ght[{index}] must increase surface to top")
    if np.any(height[0] <= _array(dataset.variables["orography"])):
        _add_error(result, f"lowest ght[{index}] must be above source orography")
    _close(result, f"two-bin dust mass[{index}]", fields["trc01"] + fields["trc02"], fields["dustq"])
    gravity = float(getattr(dataset, "gravity_m_s2", GRAVITY_M_S2))
    scale = pressure / (fields["rho"] * gravity)
    layer_height = scale * np.log(interfaces[:-1] / interfaces[1:])
    bottom_height = np.concatenate((np.zeros_like(layer_height[:1]), np.cumsum(layer_height[:-1], axis=0)), axis=0)
    expected_height = _array(dataset.variables["orography"])[None] + bottom_height
    expected_height += scale * np.log(interfaces[:-1] / pressure)
    _close(result, f"piecewise-isothermal hydrostatic ght[{index}]", height, expected_height, atol=0.03)
    column = np.sum(fields["dustq"] * (interfaces[:-1] - interfaces[1:]) / gravity, axis=0)
    _close(result, f"dust_column_mass[{index}]", fields["dust_column_mass"], column)
    ratio = fields["tau_visible610_spherical"] / np.maximum(fields["tau_visible610_mcd"], 1.0e-30)
    diagnostics = {"record": index, "dust_column_mass_min": float(np.min(column)),
        "dust_column_mass_max": float(np.max(column)), "optical_approximation_ratio_min": float(np.min(ratio)),
        "optical_approximation_ratio_max": float(np.max(ratio))}
    if "radius_clipped" in dataset.variables:
        clipped = _array(dataset.variables["radius_clipped"], index)
        if _finite(result, f"radius_clipped[{index}]", clipped):
            if np.any((clipped != 0.0) & (clipped != 1.0)):
                _add_error(result, f"radius_clipped[{index}] must be a binary mask")
            diagnostics["radius_clipped_cell_fraction"] = float(np.mean(clipped))
            if np.any(clipped > 0.0):
                fraction = diagnostics["radius_clipped_cell_fraction"]
                result["warnings"].append(f"record {index}: {fraction:.2%} of layer cells exceed bin-radius range;")
                result["warnings"][-1] += " conserved dust mass does not imply matched source effective radius"
            layer_mass = fields["dustq"] * (interfaces[:-1] - interfaces[1:]) / gravity
            mass_total = float(np.sum(layer_mass))
            clipped_mass_fraction = float(np.sum(layer_mass * clipped)) / mass_total if mass_total > 0.0 else 0.0
            diagnostics["radius_clipped_mass_fraction_unweighted_grid"] = clipped_mass_fraction
    for name in dataset.variables:
        if name.startswith("rms_"):
            values = _array(dataset.variables[name], index)
            if _finite(result, f"{name}[{index}]", values) and np.any(values < 0.0):
                _add_error(result, f"{name}[{index}] must be an unsigned nonnegative scale")
    result["sample_diagnostics"].append(diagnostics)


def _validate_file(path: Path, sample_records: int) -> tuple[dict[str, Any], np.ndarray]:
    from netCDF4 import Dataset

    result: dict[str, Any] = {"path": str(path), "errors": [], "warnings": [], "sample_diagnostics": []}
    times = np.array([], dtype=float)
    try:
        with Dataset(path, "r") as raw_dataset:
            dataset = DatasetView(raw_dataset)
            calendar_only = "pressure" not in dataset.variables and "Time" in dataset.variables
            result["calendar_only"] = calendar_only
            if calendar_only:
                required = ("Time", "Times", "solar_longitude", "mean_solar_time")
                missing = [name for name in required if name not in dataset.variables]
                if missing:
                    _add_error(result, f"calendar file lacks variables: {', '.join(missing)}")
                else:
                    times = _check_time(dataset, result)
                    result["warnings"].append("calendar-only file: no physical background was validated")
            elif _check_metadata(dataset, result):
                times = _check_time(dataset, result)
                _check_grid(dataset, result)
                indices = np.unique(np.linspace(0, times.size - 1, min(sample_records, times.size), dtype=int))
                for index in indices:
                    _check_record(dataset, result, int(index))
                result["sampled_records"] = indices.tolist()
    except (OSError, ValueError, TypeError, IndexError, ZeroDivisionError) as error:
        _add_error(result, f"validation could not complete: {error}")
    result["valid"] = not result["errors"]
    return result, times


def validate(path: str | Path, sample_records: int = 3,
    expect_full_year: bool = False) -> dict[str, Any]:
    """Validate a file or flat directory; never modify data or create logs.

    Samples include the first/last records by default. ``expect_full_year``
    requires 8028 nonoverlapping MY36 records, from sol 1 00 to sol 669 22.
    A clean subset is valid without this flag, with its coverage stated.
    """
    if isinstance(sample_records, bool) or int(sample_records) != sample_records or sample_records < 1:
        raise ValueError("sample_records must be a positive integer")
    source = Path(path).expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(source)
    files = sorted(source.glob("*.nc")) if source.is_dir() else [source]
    if not files:
        raise ValueError(f"no NetCDF files found in {source}")
    result: dict[str, Any] = {"valid": True, "path": str(source), "files": [], "errors": [], "warnings": []}
    coordinates = []
    calendar_coordinates = []
    for filename in files:
        report, times = _validate_file(filename, int(sample_records))
        result["files"].append(report)
        if report.get("calendar_only", False):
            calendar_coordinates.append(times)
        else:
            coordinates.append(times)
        result["errors"].extend(f"{filename.name}: {message}" for message in report["errors"])
        result["warnings"].extend(f"{filename.name}: {message}" for message in report["warnings"])
    active_coordinates = coordinates or calendar_coordinates
    all_times = np.sort(np.concatenate(active_coordinates)) if active_coordinates else np.array([], dtype=float)
    result["time_count"] = int(all_times.size)
    if all_times.size:
        result["first_elapsed_sol"] = float(all_times[0])
        result["last_elapsed_sol"] = float(all_times[-1])
        if all_times.size > 1 and not np.allclose(np.diff(all_times), 1.0 / SAMPLES_PER_SOL,
            rtol=0.0, atol=1.0e-7):
            _add_error(result, "combined file coverage has duplicate, overlapping, or missing two-hour records")
    if expect_full_year:
        expected = np.arange(MY36_SOL_COUNT * SAMPLES_PER_SOL) / SAMPLES_PER_SOL
        if all_times.shape != expected.shape or not np.allclose(all_times, expected, rtol=0.0, atol=1.0e-7):
            _add_error(result, "complete MY36 requires 8028 records spanning sol 1 00 to sol 669 22")
    result["valid"] = not result["errors"]
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", help="MCD background NetCDF file or directory")
    parser.add_argument("--sample-records", type=int, default=3)
    parser.add_argument("--expect-full-year", action="store_true")
    options = parser.parse_args()
    report = validate(options.path, options.sample_records, options.expect_full_year)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
