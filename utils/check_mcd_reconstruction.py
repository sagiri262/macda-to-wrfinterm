"""Measure monthly back-averaging residuals of the continuous MCD background.

Uniform two-Martian-hour samples are averaged within twelve 30-degree Ls
bins, on unchanged native hybrid-level indices. This is an interpolation and
sampling diagnostic, not an estimate of true-state error or a mean correction.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from netCDF4 import Dataset
import numpy as np

from macda2wrf.grid import orient_lat_lon
from macda2wrf.mcd_pipeline import new_dataset
from macda2wrf.mcd_reader import McdReader, read_array
from macda2wrf.mcd_time import MCD_ORBIT_SOLS, SI_SECONDS_PER_SOL, iter_mcd_times


DEFAULT_INPUT_DIR = "/home/zy/WRF/DATA/mcd_reanalysis/background"
DEFAULT_SOURCE_DIR = "/home/zy/WRF/DATA/MCD-archive/MY36"
DEFAULT_MOUNTAIN_FILE = "/mnt/c/Mars_EMM/MCD/MCD_6.2/data/mountain.nc"
DEFAULT_OUTPUT_FILE = "/home/zy/WRF/DATA/mcd_reanalysis/MY36_monthly_reconstruction_residual.nc"
FIELDS = {"temperature": ("temp", "K"), "u": ("u", "m s-1"), "v": ("v", "m s-1"),
    "ps": ("ps", "Pa"), "tau_visible610_mcd": ("tau_pref_gcm", "1"), "dustq": ("dustq", "kg kg-1")}
SAMPLE_COUNT = 8028


def check_coverage(input_dir: str | Path) -> list[dict[str, Any]]:
    """Check all MY36 dates before any atmospheric arrays are read."""
    directory = Path(input_dir).expanduser().resolve()
    if not directory.is_dir():
        raise NotADirectoryError(directory)
    files = sorted(directory.glob("*.nc"))
    if not files:
        raise ValueError(f"no background NetCDF files in {directory}")
    points = list(iter_mcd_times(36, 2))
    expected_time = np.asarray([point.elapsed_seconds for point in points])
    expected_dates = np.asarray([point.wrf_date for point in points])
    manifest = []
    all_times = []
    all_dates = []
    for path in files:
        with Dataset(path) as dataset:
            for name in ("time", "wrf_date", "ls", "orbit_wrapped"):
                if name not in dataset.variables:
                    raise ValueError(f"{path.name}: missing time metadata {name}")
            if str(getattr(dataset["time"], "units", "")) not in ("s", "seconds"):
                raise ValueError(f"{path.name}: time must be SI elapsed seconds")
            elapsed = read_array(dataset["time"])
            labels = np.asarray(dataset["wrf_date"][:], dtype=str)
            ls = read_array(dataset["ls"])
            wrapped_raw = read_array(dataset["orbit_wrapped"])
            shape = elapsed.shape
            if elapsed.ndim != 1 or elapsed.size == 0 or labels.shape != shape or ls.shape != shape:
                raise ValueError(f"{path.name}: inconsistent one-dimensional time metadata")
            if wrapped_raw.shape != shape or np.any((wrapped_raw != 0) & (wrapped_raw != 1)):
                raise ValueError(f"{path.name}: orbit_wrapped must be a time-indexed binary flag")
            if not np.all(np.isfinite(elapsed)) or np.any(np.diff(elapsed) <= 0):
                raise ValueError(f"{path.name}: time must be finite and strictly increasing")
            if not np.all(np.isfinite(ls)) or np.any(ls < 0) or np.any(ls >= 360):
                raise ValueError(f"{path.name}: ls must lie within [0, 360)")
            wrapped = wrapped_raw.astype(bool)
            expected_wrap = elapsed / SI_SECONDS_PER_SOL >= MCD_ORBIT_SOLS
            if not np.array_equal(wrapped, expected_wrap):
                raise ValueError(f"{path.name}: orbit_wrapped differs from the documented orbital boundary")
        manifest.append({"path": path, "time": elapsed, "dates": labels, "ls": ls, "wrapped": wrapped})
        all_times.append(elapsed)
        all_dates.append(labels)
    times = np.concatenate(all_times)
    dates = np.concatenate(all_dates)
    order = np.argsort(times)
    if times.size != SAMPLE_COUNT or not np.allclose(times[order], expected_time, rtol=0, atol=1.0e-5):
        raise ValueError("complete MY36 requires exactly 8028 unique two-hour samples; coverage is missing or overlaps")
    if not np.array_equal(dates[order], expected_dates):
        raise ValueError("WRF date strings do not match the complete MY36-00001 00 through MY36-00669 22 calendar")
    return sorted(manifest, key=lambda item: item["time"][0])


def _grid_matches(dataset: Any, reader: McdReader) -> None:
    for name, expected in (("lat", np.sort(reader.latitude)), ("lon", np.sort(reader.longitude))):
        actual = read_array(dataset[name])
        if actual.shape != expected.shape or not np.allclose(actual, expected, rtol=0, atol=2.0e-5):
            raise ValueError(f"{dataset.filepath()}: {name} differs from the oriented source grid")
    if len(dataset.dimensions["level"]) != reader.aps.size:
        raise ValueError(f"{dataset.filepath()}: hybrid level count differs from source")
    for name, expected in (("aps", reader.aps), ("bps", reader.bps)):
        if name in dataset.variables:
            actual = read_array(dataset[name])
            if actual.shape != expected.shape or not np.allclose(actual, expected, rtol=1.0e-7, atol=1.0e-10):
                raise ValueError(f"{dataset.filepath()}: {name} differs from source hybrid coordinates")


def _accumulate(manifest: list[dict[str, Any]], reader: McdReader) -> tuple[dict[str, np.ndarray], np.ndarray, int]:
    shape = (len(reader.latitude), len(reader.longitude))
    sums = {name: np.zeros((12,) + ((len(reader.aps),) if name not in ("ps", "tau_visible610_mcd") else ())
        + shape, dtype=np.float64) for name in FIELDS}
    counts = np.zeros(12, dtype=np.int64)
    excluded = 0
    for item in manifest:
        valid = ~item["wrapped"]
        excluded += int(np.sum(~valid))
        months = np.floor(item["ls"] / 30.0).astype(int)
        for month in range(12):
            counts[month] += int(np.sum(valid & (months == month)))
        with Dataset(item["path"]) as dataset:
            _grid_matches(dataset, reader)
            product = str(getattr(dataset, "product_type", ""))
            if product != "mcd_interpolated_background":
                raise ValueError(f"{item['path'].name}: expected unperturbed MCD background; got {product!r}")
            for name, (source_name, units) in FIELDS.items():
                if name not in dataset.variables:
                    raise ValueError(f"{item['path'].name}: missing {name}")
                variable = dataset[name]
                if str(getattr(variable, "units", "")) != units:
                    raise ValueError(f"{item['path'].name}: inconsistent units for {name}")
                values = read_array(variable)
                if values.shape != (valid.size,) + sums[name].shape[1:] or not np.all(np.isfinite(values)):
                    raise ValueError(f"{item['path'].name}: {name} has incompatible dimensions or missing values")
                for month in np.unique(months[valid]):
                    sums[name][month] += np.sum(values[valid & (months == month)], axis=0)
    if np.any(counts == 0):
        raise ValueError("one or more 30-degree Ls bins have no unwrapped background samples")
    means = {name: values / counts.reshape((12,) + (1,) * (values.ndim - 1)) for name, values in sums.items()}
    return means, counts, excluded


def check_reconstruction(input_dir: str | Path = DEFAULT_INPUT_DIR,
    source_dir: str | Path = DEFAULT_SOURCE_DIR, mountain_file: str | Path = DEFAULT_MOUNTAIN_FILE,
    output_file: str | Path = DEFAULT_OUTPUT_FILE) -> Path:
    """Publish source means, reconstructed means and reconstructed-minus-source residuals."""
    output_path = Path(output_file).expanduser().resolve()
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite {output_path}")
    manifest = check_coverage(input_dir)
    reader = McdReader(source_dir, mountain_file, year=36)
    means, counts, excluded = _accumulate(manifest, reader)
    with new_dataset(output_path) as output:
        output.product_type = "mcd_monthly_reconstruction_residual"
        output.definition = "uniform 2h Ls-bin mean minus arithmetic mean of 12 source diurnal samples"
        output.not_truth_error = "interpolation/seasonal-bin back-averaging residual; not observed truth error"
        output.mean_correction = "none: no forced restoration of source monthly means"
        output.vertical_comparison = "same native hybrid level index; not a common isobaric surface"
        output.time_weighting = "equal 2-Martian-hour samples; exclude orbit_wrapped records"
        output.source_directory = str(Path(source_dir).resolve())
        output.background_directory = str(Path(input_dir).resolve())
        output.checked_calendar_records = SAMPLE_COUNT
        output.included_records = int(np.sum(counts))
        output.excluded_wrapped_count = excluded
        output.createDimension("month", 12)
        output.createDimension("bounds", 2)
        for name, values, units in (("lat", np.sort(reader.latitude), "degrees_north"),
            ("lon", np.sort(reader.longitude), "degrees_east"), ("level", np.arange(reader.aps.size), "1")):
            output.createDimension(name, len(values))
            variable = output.createVariable(name, "f8", (name,))
            variable[:] = values
            variable.units = units
        variable = output.createVariable("month", "i4", ("month",))
        variable[:] = np.arange(1, 13)
        variable.long_name = "MCD 30-degree Ls seasonal month"
        variable = output.createVariable("ls_bounds", "f8", ("month", "bounds"))
        variable[:] = np.column_stack((np.arange(0, 360, 30), np.arange(30, 361, 30)))
        variable.units = "degree"
        variable = output.createVariable("month_sample_count", "i4", ("month",))
        variable[:] = counts
        variable.units = "1"
        variable = output.createVariable("excluded_wrapped_count", "i4")
        variable.assignValue(excluded)
        variable.units = "1"
        for name, values, units in (("aps", reader.aps, "Pa"), ("bps", reader.bps, "1")):
            variable = output.createVariable(name, "f8", ("level",))
            variable[:] = values
            variable.units = units
        variables = {}
        for name, values in means.items():
            dimensions = ("month", "level", "lat", "lon") if values.ndim == 4 else ("month", "lat", "lon")
            for suffix in ("source_monthly_mean", "reconstructed_monthly_mean", "residual"):
                variable = output.createVariable(f"{name}_{suffix}", "f4", dimensions,
                    zlib=True, complevel=1, shuffle=True, fill_value=np.nan)
                variable.units = FIELDS[name][1]
                variable.long_name = f"{name} {suffix.replace('_', ' ')}"
                variables[(name, suffix)] = variable
            variables[(name, "reconstructed_monthly_mean")][:] = values
        for month in range(1, 13):
            with Dataset(reader.month_path(month)) as source:
                reader.validate_grid(source)
                if len(source.dimensions["Time"]) != 12:
                    raise ValueError(f"source month {month} must have 12 diurnal samples")
                for name, expected in (("aps", reader.aps), ("bps", reader.bps), ("Time", reader.hours)):
                    if not np.array_equal(read_array(source[name]), expected):
                        raise ValueError(f"source month {month}: {name} differs from the common hybrid/time grid")
                for name, (source_name, units) in FIELDS.items():
                    source_values = read_array(source[source_name])
                    expected_shape = (12,) + means[name].shape[1:]
                    if source_values.shape != expected_shape or not np.all(np.isfinite(source_values)):
                        raise ValueError(f"source month {month}: invalid shape or missing values in {source_name}")
                    source_mean = np.mean(source_values, axis=0)
                    source_mean = orient_lat_lon(source_mean, reader.latitude, reader.longitude, "native")[0]
                    variables[(name, "source_monthly_mean")][month - 1] = source_mean
                    variables[(name, "residual")][month - 1] = means[name][month - 1] - source_mean
    return output_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default=DEFAULT_INPUT_DIR)
    parser.add_argument("--source-dir", default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--mountain-file", default=DEFAULT_MOUNTAIN_FILE)
    parser.add_argument("--output-file", default=DEFAULT_OUTPUT_FILE)
    options = parser.parse_args()
    path = check_reconstruction(options.input_dir, options.source_dir, options.mountain_file, options.output_file)
    print(f"monthly reconstruction residual: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
