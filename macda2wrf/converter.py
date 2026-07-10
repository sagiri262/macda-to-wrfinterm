"""MACDA v2.0 conversion driver."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from macda2wrf.config import MacdaConfig
from macda2wrf.grid import RegularLatLonGrid, horizontal_interp_2d, orient_lat_lon
from macda2wrf.macda_reader import MacdaReader
from macda2wrf.mars_time import make_hdate, parse_mars_date, validate_macda_time
from macda2wrf.vertical import interp_sigma_to_pressure
from macda2wrf.wrf_intermediate import (
    WrfIntermediateWriter,
    XLVL_SURFACE,
    validate_intermediate_file,
)


@dataclass(frozen=True)
class FieldSpec:
    src_v: str
    aim_v: str
    units: str
    kind: str
    transform: str
    required: bool
    enabled: bool
    desc: str


class MacdaConverter:
    def __init__(self, cfg: MacdaConfig):
        self.cfg = cfg
        self.grid = RegularLatLonGrid(
            lat_start=cfg.lat_start,
            lat_end=cfg.lat_end,
            nlat=cfg.nlat,
            lon_start=cfg.lon_start,
            lon_end=cfg.lon_end,
            nlon=cfg.nlon,
            radius_km=cfg.mars_radius_km,
        )
        self.fields = read_field_table(cfg.variable_table)

    def run(
        self,
        time_indices: Iterable[int] | None = None,
        max_times_override: int | None = None,
        dry_run: bool = False,
    ) -> list[Path]:
        written: list[Path] = []
        with MacdaReader(self.cfg.input_file) as reader:
            lat = reader.read("lat")
            lon = reader.read("lon")
            sigma = reader.read("lev")
            time_values = reader.read("time")
            ntimes = reader.size("time")
            indices = (
                list(time_indices)
                if time_indices is not None
                else self._time_indices(
                    ntimes, apply_config_limit=max_times_override is None
                )
            )
            if max_times_override is not None:
                if max_times_override < 1:
                    raise ValueError("max_times must be at least 1")
                indices = indices[:max_times_override]
            self._validate_indices(indices, ntimes)

            if dry_run:
                self._print_dry_run(reader, indices, time_values)
                return []

            for time_index in indices:
                print(f"[macda2w] processing time_index={time_index}")
                out_path = self._write_one_time(
                    reader,
                    time_index,
                    indices[0],
                    time_values,
                    lat,
                    lon,
                    sigma,
                )
                written.append(out_path)
        return written

    @staticmethod
    def _validate_indices(indices: list[int], ntimes: int) -> None:
        if not indices:
            raise ValueError("No MACDA time indices selected")
        invalid = [index for index in indices if index < 0 or index >= ntimes]
        if invalid:
            raise ValueError(
                f"time indices {invalid} outside available range 0..{ntimes - 1}"
            )

    def _time_indices(
        self, ntimes: int, apply_config_limit: bool = True
    ) -> list[int]:
        start = self.cfg.start_index
        end = self.cfg.end_index if self.cfg.end_index is not None else ntimes - 1
        end = min(end, ntimes - 1)
        if start < 0 or start >= ntimes:
            raise ValueError(
                f"start_index {start} outside available time range 0..{ntimes - 1}"
            )
        if end < start:
            raise ValueError(f"end_index {end} is before start_index {start}")
        indices = list(range(start, end + 1))
        if apply_config_limit and self.cfg.max_times is not None:
            indices = indices[: self.cfg.max_times]
        return indices

    def _write_one_time(
        self,
        reader: MacdaReader,
        time_index: int,
        first_time_index: int,
        time_values: np.ndarray,
        lat: np.ndarray,
        lon: np.ndarray,
        sigma: np.ndarray,
    ) -> Path:
        mars_date = reader.read_time_string(time_index)
        hdate, stamp = make_hdate(mars_date, self.cfg.hdate_strategy)
        parsed_time = parse_mars_date(mars_date)
        if self.cfg.validate_time_alignment:
            validate_macda_time(
                parsed_time,
                float(time_values[time_index]),
                tolerance_seconds=self.cfg.time_tolerance_seconds,
            )
        output_path = self.cfg.output_root / f"{self.cfg.output_prefix}:{stamp}"

        psfc_raw = reader.read("psurf", time_index=time_index)
        psfc, src_lats, src_lons = orient_lat_lon(
            psfc_raw, lat, lon, self.cfg.lon_convention
        )
        del src_lats, src_lons

        with WrfIntermediateWriter(
            output_path,
            self.grid,
            hdate=hdate,
            map_source=self.cfg.map_source,
            xfcst=self._xfcst_hours(time_values, time_index, first_time_index),
        ) as writer:
            for field in self.fields:
                if not field.enabled:
                    continue
                if not field.required and not self.cfg.emit_optional_fields:
                    continue
                self._write_field(
                    reader, writer, field, time_index, lat, lon, sigma, psfc
                )
        validate_intermediate_file(
            output_path,
            expected_hdate=hdate,
            expected_shape=(self.grid.nlat, self.grid.nlon),
            expected_xfcst=self._xfcst_hours(
                time_values, time_index, first_time_index
            ),
            required_fields={
                field.aim_v
                for field in self.fields
                if field.enabled and field.required
            },
        )
        return output_path

    def _write_field(
        self,
        reader: MacdaReader,
        writer: WrfIntermediateWriter,
        field: FieldSpec,
        time_index: int,
        lat: np.ndarray,
        lon: np.ndarray,
        sigma: np.ndarray,
        psfc: np.ndarray,
    ) -> None:
        if field.kind == "3d_pressure":
            for plev in self.cfg.plev_pa:
                slab = np.full(
                    (self.grid.nlat, self.grid.nlon), plev, dtype=np.float32
                )
                writer.write_field(
                    field.aim_v, slab, field.units, field.desc, xlvl=plev
                )
            return

        if field.kind == "3d_constant":
            for plev in self.cfg.plev_pa:
                slab = np.zeros((self.grid.nlat, self.grid.nlon), dtype=np.float32)
                writer.write_field(field.aim_v, slab, field.units, field.desc, xlvl=plev)
            return

        if not reader.has_var(field.src_v):
            if field.required:
                raise KeyError(f"Required MACDA variable missing: {field.src_v}")
            print(f"[macda2w] skip optional missing variable {field.src_v}")
            return

        if field.kind == "2d":
            raw = reader.read(field.src_v, time_index=time_index)
            slab, src_lats, src_lons = orient_lat_lon(
                raw, lat, lon, self.cfg.lon_convention
            )
            slab = self._transform_2d(field, slab, psfc)
            out = horizontal_interp_2d(slab, src_lats, src_lons, self.grid)
            self._ensure_finite(field, out)
            writer.write_field(
                field.aim_v,
                out,
                field.units,
                field.desc,
                xlvl=XLVL_SURFACE,
                is_wind_earth_rel=int(field.aim_v in {"UU", "VV"}),
            )
            return

        if field.kind == "3d_sigma":
            raw = reader.read(field.src_v, time_index=time_index)
            volume, src_lats, src_lons = orient_lat_lon(
                raw, lat, lon, self.cfg.lon_convention
            )
            volume = self._transform_3d(field, volume)
            pvol = interp_sigma_to_pressure(volume, psfc, sigma, self.cfg.plev_pa)
            for idx, plev in enumerate(self.cfg.plev_pa):
                out = horizontal_interp_2d(pvol[idx], src_lats, src_lons, self.grid)
                self._ensure_finite(field, out)
                writer.write_field(
                    field.aim_v,
                    out,
                    field.units,
                    field.desc,
                    xlvl=plev,
                    is_wind_earth_rel=int(field.aim_v in {"UU", "VV"}),
                )
            return

        raise ValueError(f"Unsupported field kind {field.kind} for {field.aim_v}")

    @staticmethod
    def _ensure_finite(field: FieldSpec, data: np.ndarray) -> None:
        bad = int(np.size(data) - np.isfinite(data).sum())
        if bad:
            raise ValueError(
                f"Field {field.aim_v} contains {bad} non-finite values after conversion"
            )

    def _transform_2d(
        self, field: FieldSpec, slab: np.ndarray, psfc: np.ndarray
    ) -> np.ndarray:
        if field.transform == "identity":
            return slab
        if field.transform == "dust_to_tau_7mb":
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.where(
                    psfc > 0.0,
                    slab / psfc * self.cfg.tau_reference_pressure_pa,
                    np.nan,
                )
        raise ValueError(f"Unsupported 2-D transform {field.transform}")

    def _transform_3d(self, field: FieldSpec, volume: np.ndarray) -> np.ndarray:
        if field.transform == "identity":
            return volume
        if field.transform == "geop_to_height":
            return volume / self.cfg.mars_gravity
        raise ValueError(f"Unsupported 3-D transform {field.transform}")

    def _xfcst_hours(
        self,
        time_values: np.ndarray,
        time_index: int,
        first_time_index: int,
    ) -> float:
        if not self.cfg.xfcst_from_start:
            return 0.0
        return float((time_values[time_index] - time_values[first_time_index]) * 24.0)

    def _print_dry_run(
        self,
        reader: MacdaReader,
        indices: list[int],
        time_values: np.ndarray,
    ) -> None:
        print(f"input_file={self.cfg.input_file}")
        print(f"reader_backend={reader.backend}")
        print(f"output_root={self.cfg.output_root}")
        print(f"time_indices={indices}")
        print(f"target_grid={self.grid.nlat}x{self.grid.nlon}")
        print(f"plev_pa={self.cfg.plev_pa}")
        for index in dict.fromkeys((indices[0], indices[-1])):
            mars_date = reader.read_time_string(index)
            parsed = parse_mars_date(mars_date)
            if self.cfg.validate_time_alignment:
                validate_macda_time(
                    parsed,
                    float(time_values[index]),
                    tolerance_seconds=self.cfg.time_tolerance_seconds,
                )
            print(
                f"time[{index}]={float(time_values[index]):.12g} sol "
                f"{mars_date} -> {parsed.wrf_date}"
            )
        print("fields:")
        for field in self.fields:
            if field.enabled and (field.required or self.cfg.emit_optional_fields):
                status = "required" if field.required else "optional"
                print(
                    f"  {field.aim_v:10s} <- {field.src_v:12s} "
                    f"{field.kind:13s} {status}"
                )


def read_field_table(path: str | Path) -> list[FieldSpec]:
    fields: list[FieldSpec] = []
    with Path(path).open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            field = FieldSpec(
                src_v=row["src_v"].strip(),
                aim_v=row["aim_v"].strip(),
                units=row["units"].strip(),
                kind=row["kind"].strip(),
                transform=row["transform"].strip(),
                required=_csv_bool(row["required"]),
                enabled=_csv_bool(row["enabled"]),
                desc=row["desc"].strip(),
            )
            if not field.src_v or not field.aim_v:
                raise ValueError(f"Empty src_v/aim_v in variable table {path}")
            if len(field.aim_v) > 9:
                raise ValueError(
                    f"WRF intermediate field name exceeds 9 characters: {field.aim_v}"
                )
            if field.kind not in {
                "2d",
                "3d_sigma",
                "3d_pressure",
                "3d_constant",
            }:
                raise ValueError(f"Unsupported kind {field.kind!r} for {field.aim_v}")
            fields.append(field)
    if not fields:
        raise ValueError(f"Variable table contains no fields: {path}")
    return fields


def _csv_bool(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "on"}
