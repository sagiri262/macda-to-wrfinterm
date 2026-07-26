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
from macda2wrf.mars_time import (ANCHORED_HDATE_STRATEGIES, LabelledRecord,
                                 MacdaRecord, MARSWRF_FIXED_YEAR_SOLS,
                                 build_label_sequence, short_mars_years_left,
                                 validate_macda_time)
from macda2wrf.vertical import interp_sigma_to_pressure
from macda2wrf.wrf_intermediate import (WrfIntermediateWriter, XLVL_SURFACE,
                                        validate_intermediate_file)


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
        self.grid = RegularLatLonGrid(lat_start=cfg.lat_start, lat_end=cfg.lat_end,
                                      nlat=cfg.nlat, lon_start=cfg.lon_start,
                                      lon_end=cfg.lon_end, nlon=cfg.nlon,
                                      radius_km=cfg.mars_radius_km)
        self.fields = read_field_table(cfg.variable_table)

    def run(self, time_indices: Iterable[int] | None = None,
            max_times_override: int | None = None,
            dry_run: bool = False) -> list[Path]:
        written: list[Path] = []
        with MacdaReader(self.cfg.input_file) as reader:
            lat = reader.read("lat")
            lon = reader.read("lon")
            sigma = reader.read("lev")
            time_values = reader.read("time")
            ntimes = reader.size("time")
            anchor_index = self._anchor_index(ntimes)
            indices = self._select_indices(ntimes, anchor_index, time_indices,
                                           max_times_override)
            labelled = self._build_labels(reader, indices, time_values, anchor_index)

            if dry_run:
                self._print_dry_run(reader, labelled, time_values, anchor_index)
                return []

            for record in labelled:
                print(f"[macda2w] processing time_index={record.time_index}")
                out_path = self._write_one_time(reader, record, time_values,
                                                anchor_index, lat, lon, sigma)
                written.append(out_path)
        return written

    def _anchor_index(self, ntimes: int) -> int:
        """Return the stable label and XFCST anchor, which is ``start_index``.

        The anchor never depends on which time indices a single command selects,
        so a partial rerun writes exactly the filenames a full run writes.
        """

        anchor_index = self.cfg.start_index
        if anchor_index < 0 or anchor_index >= ntimes:
            raise ValueError(f"start_index {anchor_index} outside available time range "
                             f"0..{ntimes - 1}; it is the label and XFCST anchor")
        return anchor_index

    def _select_indices(self, ntimes: int, anchor_index: int,
                        time_indices: Iterable[int] | None,
                        max_times_override: int | None) -> list[int]:
        if time_indices is None:
            apply_config_limit = max_times_override is None
            indices = self._time_indices(ntimes, apply_config_limit=apply_config_limit)
        else:
            requested = [int(index) for index in time_indices]
            indices = sorted(dict.fromkeys(requested))
            if indices != requested:
                print(f"[macda2w] explicit time indices sorted and deduplicated: "
                      f"{indices}")
        if max_times_override is not None:
            if max_times_override < 1:
                raise ValueError("max_times must be at least 1")
            indices = indices[:max_times_override]
        self._validate_indices(indices, ntimes, anchor_index)
        return indices

    @staticmethod
    def _validate_indices(indices: list[int], ntimes: int, anchor_index: int) -> None:
        if not indices:
            raise ValueError("No MACDA time indices selected")
        invalid = [index for index in indices if index < 0 or index >= ntimes]
        if invalid:
            raise ValueError(f"time indices {invalid} outside available range "
                             f"0..{ntimes - 1}")
        early = [index for index in indices if index < anchor_index]
        if early:
            raise ValueError(f"time indices {early} come before the anchor record "
                             f"start_index={anchor_index}, so XFCST would be negative; "
                             f"lower start_index in the config to move the anchor")

    def _time_indices(self, ntimes: int, apply_config_limit: bool = True) -> list[int]:
        start = self.cfg.start_index
        end = self.cfg.end_index if self.cfg.end_index is not None else ntimes - 1
        end = min(end, ntimes - 1)
        if start < 0 or start >= ntimes:
            raise ValueError(f"start_index {start} outside available time range "
                             f"0..{ntimes - 1}")
        if end < start:
            raise ValueError(f"end_index {end} is before start_index {start}")
        indices = list(range(start, end + 1))
        if apply_config_limit and self.cfg.max_times is not None:
            indices = indices[: self.cfg.max_times]
        return indices

    def _build_labels(self, reader: MacdaReader, indices: list[int],
                      time_values: np.ndarray,
                      anchor_index: int) -> list[LabelledRecord]:
        """Label the whole selection against the anchor and validate the sequence."""

        anchor_date = reader.read_time_string(anchor_index)
        records = [MacdaRecord(time_index=index,
                               mars_date=reader.read_time_string(index),
                               time_sols=float(time_values[index]))
                   for index in indices]
        tolerance = self.cfg.time_tolerance_seconds
        labelled = build_label_sequence(records, strategy=self.cfg.hdate_strategy,
                                        anchor_date=anchor_date,
                                        tolerance_seconds=tolerance)
        short_years = short_mars_years_left(labelled)
        if short_years and self.cfg.hdate_strategy not in ANCHORED_HDATE_STRATEGIES:
            print(f"[macda2w] warning: this selection leaves MY {short_years}, which is "
                  f"shorter than the fixed {MARSWRF_FIXED_YEAR_SOLS}-sol MarsWRF year; "
                  f"metgrid will ask for a sol these labels do not contain")
        return labelled

    def _write_one_time(self, reader: MacdaReader, record: LabelledRecord,
                        time_values: np.ndarray, anchor_index: int, lat: np.ndarray,
                        lon: np.ndarray, sigma: np.ndarray) -> Path:
        time_index = record.time_index
        hdate = record.hdate
        stamp = record.filename_stamp
        if record.relabelled:
            print(f"[macda2w] {self.cfg.hdate_strategy} relabel "
                  f"{record.mars_date.filename_stamp} -> {stamp}")
        if self.cfg.validate_time_alignment:
            validate_macda_time(record.mars_date, float(time_values[time_index]),
                                tolerance_seconds=self.cfg.time_tolerance_seconds)
        output_path = self.cfg.output_root / f"{self.cfg.output_prefix}:{stamp}"
        xfcst = self._xfcst_hours(time_values, time_index, anchor_index)

        psfc_raw = reader.read("psurf", time_index=time_index)
        psfc, src_lats, src_lons = orient_lat_lon(psfc_raw, lat, lon,
                                                  self.cfg.lon_convention)
        del src_lats, src_lons

        with WrfIntermediateWriter(output_path, self.grid, hdate=hdate,
                                   map_source=self.cfg.map_source,
                                   xfcst=xfcst) as writer:
            for field in self.fields:
                if not field.enabled:
                    continue
                if not field.required and not self.cfg.emit_optional_fields:
                    continue
                self._write_field(reader, writer, field, time_index, lat, lon, sigma,
                                  psfc)
        required_fields = {field.aim_v for field in self.fields
                           if field.enabled and field.required}
        validate_intermediate_file(output_path, expected_hdate=hdate,
                                   expected_shape=(self.grid.nlat, self.grid.nlon),
                                   expected_xfcst=xfcst,
                                   required_fields=required_fields)
        return output_path

    def _write_field(self, reader: MacdaReader, writer: WrfIntermediateWriter,
                     field: FieldSpec, time_index: int, lat: np.ndarray,
                     lon: np.ndarray, sigma: np.ndarray, psfc: np.ndarray) -> None:
        if field.kind == "3d_pressure":
            for plev in self.cfg.plev_pa:
                slab = np.full((self.grid.nlat, self.grid.nlon), plev,
                               dtype=np.float32)
                writer.write_field(field.aim_v, slab, field.units, field.desc,
                                   xlvl=plev)
            return

        if field.kind == "3d_constant":
            for plev in self.cfg.plev_pa:
                slab = np.zeros((self.grid.nlat, self.grid.nlon), dtype=np.float32)
                writer.write_field(field.aim_v, slab, field.units, field.desc,
                                   xlvl=plev)
            return

        if not reader.has_var(field.src_v):
            if field.required:
                raise KeyError(f"Required MACDA variable missing: {field.src_v}")
            print(f"[macda2w] skip optional missing variable {field.src_v}")
            return

        if field.kind == "2d":
            raw = reader.read(field.src_v, time_index=time_index)
            slab, src_lats, src_lons = orient_lat_lon(raw, lat, lon,
                                                      self.cfg.lon_convention)
            slab = self._transform_2d(field, slab, psfc)
            out = horizontal_interp_2d(slab, src_lats, src_lons, self.grid)
            self._ensure_finite(field, out)
            writer.write_field(field.aim_v, out, field.units, field.desc,
                               xlvl=XLVL_SURFACE,
                               is_wind_earth_rel=int(field.aim_v in {"UU", "VV"}))
            return

        if field.kind == "3d_sigma":
            raw = reader.read(field.src_v, time_index=time_index)
            volume, src_lats, src_lons = orient_lat_lon(raw, lat, lon,
                                                        self.cfg.lon_convention)
            volume = self._transform_3d(field, volume)
            pvol = interp_sigma_to_pressure(volume, psfc, sigma, self.cfg.plev_pa)
            for idx, plev in enumerate(self.cfg.plev_pa):
                out = horizontal_interp_2d(pvol[idx], src_lats, src_lons, self.grid)
                self._ensure_finite(field, out)
                writer.write_field(field.aim_v, out, field.units, field.desc,
                                   xlvl=plev,
                                   is_wind_earth_rel=int(field.aim_v in {"UU", "VV"}))
            return

        raise ValueError(f"Unsupported field kind {field.kind} for {field.aim_v}")

    @staticmethod
    def _ensure_finite(field: FieldSpec, data: np.ndarray) -> None:
        bad = int(np.size(data) - np.isfinite(data).sum())
        if bad:
            raise ValueError(f"Field {field.aim_v} contains {bad} non-finite values "
                             f"after conversion")

    def _transform_2d(self, field: FieldSpec, slab: np.ndarray,
                      psfc: np.ndarray) -> np.ndarray:
        if field.transform == "identity":
            return slab
        if field.transform == "dust_to_tau_7mb":
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.where(psfc > 0.0,
                                slab / psfc * self.cfg.tau_reference_pressure_pa,
                                np.nan)
        raise ValueError(f"Unsupported 2-D transform {field.transform}")

    def _transform_3d(self, field: FieldSpec, volume: np.ndarray) -> np.ndarray:
        if field.transform == "identity":
            return volume
        if field.transform == "geop_to_height":
            return volume / self.cfg.mars_gravity
        raise ValueError(f"Unsupported 3-D transform {field.transform}")

    def _xfcst_hours(self, time_values: np.ndarray, time_index: int,
                     anchor_index: int) -> float:
        """Forecast hours since the anchor record, in Martian hours."""

        if not self.cfg.xfcst_from_start:
            return 0.0
        hours = float((time_values[time_index] - time_values[anchor_index]) * 24.0)
        if hours < 0.0:
            raise ValueError(f"XFCST for time index {time_index} is {hours:.6g} h, "
                             f"before the anchor record start_index={anchor_index}")
        return hours

    def _print_dry_run(self, reader: MacdaReader, labelled: list[LabelledRecord],
                       time_values: np.ndarray, anchor_index: int) -> None:
        anchor_date = reader.read_time_string(anchor_index)
        relabelled = [record for record in labelled if record.relabelled]
        print(f"input_file={self.cfg.input_file}")
        print(f"reader_backend={reader.backend}")
        print(f"output_root={self.cfg.output_root}")
        print(f"target_grid={self.grid.nlat}x{self.grid.nlon}")
        print(f"plev_pa={self.cfg.plev_pa}")
        print(f"hdate_strategy={self.cfg.hdate_strategy}")
        print(f"label_anchor=time[{anchor_index}] {anchor_date}")
        print(f"records={len(labelled)} relabelled={len(relabelled)} "
              f"first={labelled[0].hdate[:19]} last={labelled[-1].hdate[:19]}")
        print("times:")
        for record in labelled:
            if self.cfg.validate_time_alignment:
                validate_macda_time(record.mars_date, record.time_sols,
                                    tolerance_seconds=self.cfg.time_tolerance_seconds)
            xfcst = self._xfcst_hours(time_values, record.time_index, anchor_index)
            mark = " relabelled" if record.relabelled else ""
            print(f"  time[{record.time_index}]={record.time_sols:.12g} sol "
                  f"{record.mars_date.macda_string} -> {record.mars_date.wrf_date} "
                  f"-> label {record.hdate[:19]} xfcst={xfcst:.6g} h{mark}")
        print("fields:")
        for field in self.fields:
            if field.enabled and (field.required or self.cfg.emit_optional_fields):
                status = "required" if field.required else "optional"
                print(f"  {field.aim_v:10s} <- {field.src_v:12s} "
                      f"{field.kind:13s} {status}")


def read_field_table(path: str | Path) -> list[FieldSpec]:
    fields: list[FieldSpec] = []
    with Path(path).open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            field = FieldSpec(src_v=row["src_v"].strip(), aim_v=row["aim_v"].strip(),
                              units=row["units"].strip(), kind=row["kind"].strip(),
                              transform=row["transform"].strip(),
                              required=_csv_bool(row["required"]),
                              enabled=_csv_bool(row["enabled"]),
                              desc=row["desc"].strip())
            if not field.src_v or not field.aim_v:
                raise ValueError(f"Empty src_v/aim_v in variable table {path}")
            if len(field.aim_v) > 9:
                raise ValueError(f"WRF intermediate field name exceeds 9 characters: "
                                 f"{field.aim_v}")
            if field.kind not in {"2d", "3d_sigma", "3d_pressure", "3d_constant"}:
                raise ValueError(f"Unsupported kind {field.kind!r} for {field.aim_v}")
            fields.append(field)
    if not fields:
        raise ValueError(f"Variable table contains no fields: {path}")
    return fields


def _csv_bool(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "on"}
