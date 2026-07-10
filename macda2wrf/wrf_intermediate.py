"""WRF intermediate format writer."""

from __future__ import annotations

from pathlib import Path
import struct
from dataclasses import dataclass

import numpy as np
from scipy.io import FortranEOFError, FortranFile

from macda2wrf.grid import RegularLatLonGrid


XLVL_SURFACE = 200100.0


class WrfIntermediateWriter:
    def __init__(
        self,
        path: str | Path,
        grid: RegularLatLonGrid,
        hdate: str,
        map_source: str,
        xfcst: float = 0.0,
    ):
        self.path = Path(path)
        self.grid = grid
        self.hdate = hdate
        self.map_source = map_source
        self.xfcst = float(xfcst)
        self._fh = None

    def __enter__(self) -> "WrfIntermediateWriter":
        if len(self.hdate) != 24:
            raise ValueError(
                f"WRF intermediate HDATE must be 24 characters: {self.hdate!r}"
            )
        self.hdate.encode("ascii")
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = FortranFile(self.path, "w", header_dtype=np.dtype(">u4"))
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if self._fh is not None:
            self._fh.close()
        self._fh = None

    def write_field(
        self,
        field: str,
        slab: np.ndarray,
        units: str,
        desc: str,
        xlvl: float = XLVL_SURFACE,
        is_wind_earth_rel: int = 0,
    ) -> None:
        if self._fh is None:
            raise RuntimeError("Writer is not open")
        slab = np.asarray(slab, dtype=np.float32)
        if slab.shape != (self.grid.nlat, self.grid.nlon):
            raise ValueError(
                f"{field} slab shape {slab.shape} does not match "
                f"target grid {(self.grid.nlat, self.grid.nlon)}"
            )

        self._fh.write_record(struct.pack(">I", 5))
        header = struct.pack(
            ">24sf32s9s25s46sfIII",
            _fixed(self.hdate, 24),
            self.xfcst,
            _fixed(self.map_source, 32),
            _fixed(field, 9),
            _fixed(units, 25),
            _fixed(desc, 46),
            float(xlvl),
            int(self.grid.nlon),
            int(self.grid.nlat),
            0,
        )
        self._fh.write_record(header)

        loc = struct.pack(
            ">8sfffff",
            _fixed("SWCORNER", 8),
            float(self.grid.lat_start),
            float(self.grid.lon_start),
            float(self.grid.deltlat),
            float(self.grid.deltlon),
            float(self.grid.radius_km),
        )
        self._fh.write_record(loc)
        self._fh.write_record(struct.pack(">I", int(is_wind_earth_rel)))
        self._fh.write_record(np.asarray(slab, dtype=">f4"))


def _fixed(value: str, width: int) -> bytes:
    return str(value).encode("ascii", errors="replace")[:width].ljust(width)


@dataclass(frozen=True)
class IntermediateRecord:
    hdate: str
    xfcst: float
    field: str
    units: str
    xlvl: float
    nx: int
    ny: int
    is_wind_earth_rel: int
    slab: np.ndarray


def read_intermediate_file(path: str | Path) -> list[IntermediateRecord]:
    """Read records written by :class:`WrfIntermediateWriter`."""

    records: list[IntermediateRecord] = []
    with FortranFile(path, "r", header_dtype=np.dtype(">u4")) as fh:
        while True:
            try:
                ifv = int(fh.read_record(dtype=">u4")[0])
            except FortranEOFError:
                break
            if ifv != 5:
                raise ValueError(f"Unsupported WRF intermediate format version {ifv}")
            header_bytes = fh.read_record(dtype="u1").tobytes()
            if len(header_bytes) != struct.calcsize(">24sf32s9s25s46sfIII"):
                raise ValueError("Invalid WRF intermediate metadata record size")
            unpacked = struct.unpack(">24sf32s9s25s46sfIII", header_bytes)
            hdate_raw, xfcst, _, field_raw, units_raw, _, xlvl, nx, ny, iproj = unpacked
            if iproj != 0:
                raise ValueError(f"Only regular lat-lon IPROJ=0 is supported, got {iproj}")
            loc_bytes = fh.read_record(dtype="u1").tobytes()
            if len(loc_bytes) != struct.calcsize(">8sfffff"):
                raise ValueError("Invalid WRF intermediate projection record size")
            wind_flag = int(fh.read_record(dtype=">u4")[0])
            slab = fh.read_record(dtype=">f4")
            if slab.size != nx * ny:
                raise ValueError(
                    f"Field {_text(field_raw)} has {slab.size} values, expected {nx * ny}"
                )
            records.append(
                IntermediateRecord(
                    hdate=_text(hdate_raw),
                    xfcst=float(xfcst),
                    field=_text(field_raw),
                    units=_text(units_raw),
                    xlvl=float(xlvl),
                    nx=int(nx),
                    ny=int(ny),
                    is_wind_earth_rel=wind_flag,
                    slab=slab.reshape((ny, nx)),
                )
            )
    return records


def validate_intermediate_file(
    path: str | Path,
    expected_hdate: str,
    expected_shape: tuple[int, int],
    required_fields: set[str],
    expected_xfcst: float | None = None,
) -> None:
    records = read_intermediate_file(path)
    if not records:
        raise ValueError(f"WRF intermediate file contains no records: {path}")
    fields = {record.field for record in records}
    missing = required_fields - fields
    if missing:
        raise ValueError(f"WRF intermediate file is missing fields: {sorted(missing)}")
    for record in records:
        if record.hdate != expected_hdate:
            raise ValueError(
                f"Field {record.field} HDATE {record.hdate!r} != {expected_hdate!r}"
            )
        if expected_xfcst is not None and not np.isclose(
            record.xfcst, expected_xfcst, rtol=0.0, atol=1.0e-4
        ):
            raise ValueError(
                f"Field {record.field} XFCST {record.xfcst} != {expected_xfcst}"
            )
        if record.slab.shape != expected_shape:
            raise ValueError(
                f"Field {record.field} shape {record.slab.shape} != {expected_shape}"
            )
        if not np.isfinite(record.slab).all():
            raise ValueError(f"Field {record.field} contains non-finite values")
        if record.field in {"UU", "VV"} and record.is_wind_earth_rel != 1:
            raise ValueError(f"Field {record.field} is not marked earth-relative")


def _text(value: bytes) -> str:
    return value.decode("ascii").rstrip(" \x00")
