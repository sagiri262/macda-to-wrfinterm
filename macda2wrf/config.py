"""Configuration handling for the MACDA to WRF intermediate converter."""

from __future__ import annotations

from configparser import ConfigParser
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


def _split_floats(value: str) -> list[float]:
    return [float(item.strip()) for item in value.split(",") if item.strip()]


def _bool(value: str, default: bool = False) -> bool:
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "on"}


def _optional_int(value: str, default: Optional[int] = None) -> Optional[int]:
    if value is None or value.strip() == "":
        return default
    return int(value)


def _resolve(path_value: str, base_dir: Path) -> Path:
    path = Path(path_value).expanduser()
    if path.is_absolute():
        return path
    return (base_dir / path).resolve()


@dataclass(frozen=True)
class MacdaConfig:
    config_path: Path
    project_root: Path
    input_file: Path
    variable_table: Path
    output_root: Path
    output_prefix: str
    start_index: int
    end_index: Optional[int]
    max_times: Optional[int]
    map_source: str
    mars_radius_km: float
    mars_gravity: float
    tau_reference_pressure_pa: float
    plev_pa: list[float]
    lat_start: float
    lat_end: float
    nlat: int
    lon_start: float
    lon_end: float
    nlon: int
    lon_convention: str
    hdate_strategy: str
    validate_time_alignment: bool
    time_tolerance_seconds: float
    emit_optional_fields: bool
    xfcst_from_start: bool


def load_config(config_path: str | Path) -> MacdaConfig:
    cfg_path = Path(config_path).expanduser().resolve()
    parser = ConfigParser()
    read_files = parser.read(cfg_path)
    if not read_files:
        raise FileNotFoundError(f"Cannot read config file: {cfg_path}")

    project_root = cfg_path.parent.parent
    input_section = parser["INPUT"]
    output_section = parser["OUTPUT"]
    mars_section = parser["MARS"] if parser.has_section("MARS") else {}

    table_value = output_section.get("variable_table", "db/MACDA-v2_CORE.csv")
    output_root_value = output_section.get("output_root", "output")

    config = MacdaConfig(
        config_path=cfg_path,
        project_root=project_root,
        input_file=_resolve(input_section["macda_file"], project_root),
        variable_table=_resolve(table_value, project_root),
        output_root=_resolve(output_root_value, project_root),
        output_prefix=output_section.get("output_prefix", "MACDA"),
        start_index=input_section.getint("start_index", fallback=0),
        end_index=_optional_int(input_section.get("end_index", fallback="")),
        max_times=_optional_int(input_section.get("max_times", fallback="")),
        map_source=output_section.get("map_source", "MACDA"),
        mars_radius_km=output_section.getfloat("mars_radius_km", fallback=3389.92),
        mars_gravity=float(mars_section.get("gravity", 3.72)),
        tau_reference_pressure_pa=float(
            mars_section.get("tau_reference_pressure_pa", 700.0)
        ),
        plev_pa=_split_floats(output_section["plev_pa"]),
        lat_start=output_section.getfloat("lat_start", fallback=-87.5),
        lat_end=output_section.getfloat("lat_end", fallback=87.5),
        nlat=output_section.getint("nlat", fallback=36),
        lon_start=output_section.getfloat("lon_start", fallback=0.0),
        lon_end=output_section.getfloat("lon_end", fallback=355.0),
        nlon=output_section.getint("nlon", fallback=72),
        lon_convention=output_section.get("lon_convention", "0_360"),
        hdate_strategy=output_section.get("hdate_strategy", "marswrf_sol"),
        validate_time_alignment=_bool(
            mars_section.get("validate_time_alignment", "true"), default=True
        ),
        time_tolerance_seconds=float(
            mars_section.get("time_tolerance_seconds", 1.0)
        ),
        emit_optional_fields=_bool(
            output_section.get("emit_optional_fields", "true"), default=True
        ),
        xfcst_from_start=_bool(
            output_section.get("xfcst_from_start", "true"), default=True
        ),
    )
    _validate_config(config)
    return config


def _validate_config(config: MacdaConfig) -> None:
    if not config.input_file.is_file():
        raise FileNotFoundError(f"MACDA input file does not exist: {config.input_file}")
    if not config.variable_table.is_file():
        raise FileNotFoundError(
            f"MACDA variable table does not exist: {config.variable_table}"
        )
    if config.nlat < 2 or config.nlon < 2:
        raise ValueError("Target nlat and nlon must both be at least 2")
    if config.lat_start >= config.lat_end or config.lon_start >= config.lon_end:
        raise ValueError("Target grid start coordinates must be below end coordinates")
    if not config.plev_pa or any(level <= 0.0 for level in config.plev_pa):
        raise ValueError("plev_pa must contain positive pressure levels in Pa")
    if len(set(config.plev_pa)) != len(config.plev_pa):
        raise ValueError("plev_pa contains duplicate pressure levels")
    if config.time_tolerance_seconds < 0.0:
        raise ValueError("time_tolerance_seconds cannot be negative")
