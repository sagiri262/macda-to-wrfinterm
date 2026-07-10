"""Small NetCDF4/HDF5 reader for MACDA files.

The MACDA v2.0 files are NetCDF4/HDF5. This module intentionally avoids
xarray-specific assumptions and uses h5py or netCDF4 when either is installed.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np


class MacdaReader:
    def __init__(self, path: str | Path):
        self.path = Path(path).expanduser().resolve()
        self._backend = None
        self._handle = None

    def __enter__(self) -> "MacdaReader":
        self.open()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def open(self) -> None:
        if self._handle is not None:
            return
        try:
            import h5py  # type: ignore

            self._handle = h5py.File(self.path, "r")
            self._backend = "h5py"
            return
        except ImportError:
            pass

        try:
            from netCDF4 import Dataset  # type: ignore

            self._handle = Dataset(self.path, "r")
            self._backend = "netCDF4"
            return
        except ImportError as exc:
            raise RuntimeError(
                "MACDA v2.0 is NetCDF4/HDF5. Install h5py or netCDF4 in this "
                "Python environment before running conversion."
            ) from exc

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()
        self._handle = None
        self._backend = None

    @property
    def backend(self) -> str:
        if self._backend is None:
            raise RuntimeError("Reader is not open")
        return self._backend

    def has_var(self, name: str) -> bool:
        self._require_open()
        if self.backend == "h5py":
            return name in self._handle
        return name in self._handle.variables

    def size(self, dim_name: str) -> int:
        arr = self.read(dim_name)
        return int(arr.shape[0])

    def read(self, name: str, time_index: int | None = None) -> np.ndarray:
        self._require_open()
        obj = self._get_obj(name)
        if time_index is None:
            data = obj[()]
        else:
            data = obj[time_index]
        if np.ma.isMaskedArray(data):
            return np.asarray(data.filled(np.nan))

        array = np.asarray(data)
        if self.backend != "h5py" or array.dtype.kind not in {"i", "u", "f"}:
            return array

        fill_value = obj.attrs.get("_FillValue")
        missing_value = obj.attrs.get("missing_value")
        scale_factor = obj.attrs.get("scale_factor")
        add_offset = obj.attrs.get("add_offset")
        metadata = (fill_value, missing_value, scale_factor, add_offset)
        if any(value is not None for value in metadata):
            raw = array
            array = raw.astype(np.float64)
            invalid = np.zeros(raw.shape, dtype=bool)
            for marker in (fill_value, missing_value):
                if marker is not None:
                    invalid |= raw == np.asarray(marker).reshape(()).item()
            if scale_factor is not None:
                array *= float(np.asarray(scale_factor).reshape(()).item())
            if add_offset is not None:
                array += float(np.asarray(add_offset).reshape(()).item())
            array[invalid] = np.nan
        return array

    def read_time_string(self, time_index: int) -> str:
        raw = self.read("Mars_date", time_index=time_index)
        return _decode_scalar(raw)

    def attrs(self, name: str) -> dict[str, Any]:
        self._require_open()
        obj = self._get_obj(name)
        if self.backend == "h5py":
            return {str(key): value for key, value in obj.attrs.items()}
        return {str(key): getattr(obj, key) for key in obj.ncattrs()}

    def _get_obj(self, name: str):
        if not self.has_var(name):
            raise KeyError(f"Variable not found in MACDA file: {name}")
        if self.backend == "h5py":
            return self._handle[name]
        return self._handle.variables[name]

    def _require_open(self) -> None:
        if self._handle is None:
            raise RuntimeError("Reader is not open")


def _decode_scalar(value: np.ndarray | bytes | str) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.ndarray):
        if value.shape == ():
            return _decode_scalar(value.item())
        if value.dtype.kind in {"S", "U"}:
            return "".join(_decode_scalar(item) for item in value.ravel())
        return str(value.item())
    return str(value)
