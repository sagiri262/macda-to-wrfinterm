#!/usr/bin/env python3
"""Build a 1-degree WPS geogrid data set from the MarsWRF Data directory."""

from __future__ import annotations

import argparse
import hashlib
import json
import struct
from pathlib import Path

import numpy as np


NX = 360
NY = 180
BORDER = 3
TILE_NAME = "00001-00360.00001-00180"
TES_WRF_LON0 = -179.229  # 1-degree mean center after the IAU 1994 -> 2000 offset.


def require_size(path: Path, expected: int) -> None:
    actual = path.stat().st_size
    if actual != expected:
        raise ValueError(f"{path}: expected {expected} bytes, found {actual}")


def check_fortran_record(path: Path, payload_bytes: int) -> None:
    require_size(path, payload_bytes + 8)
    with path.open("rb") as stream:
        leading = struct.unpack(">i", stream.read(4))[0]
        stream.seek(-4, 2)
        trailing = struct.unpack(">i", stream.read(4))[0]
    if leading != payload_bytes or trailing != payload_bytes:
        raise ValueError(
            f"{path}: invalid big-endian Fortran record markers "
            f"({leading}, {trailing}); expected {payload_bytes}"
        )


def mean_blocks(source: np.ndarray, y_factor: int, x_factor: int) -> np.ndarray:
    """Average a regular source grid without loading the whole source into RAM."""
    if source.shape != (NY * y_factor, NX * x_factor):
        raise ValueError(f"unexpected source shape {source.shape}")
    result = np.empty((NY, NX), dtype=np.float64)
    for j in range(NY):
        rows = np.asarray(source[j * y_factor : (j + 1) * y_factor], dtype=np.float64)
        result[j] = rows.reshape(y_factor, NX, x_factor).mean(axis=(0, 2))
    return result


def mean_tes_map(source: np.ndarray) -> np.ndarray:
    """Apply the MarsWRF TES longitude transform and aggregate 20 px/degree data."""
    result = np.empty((NY, NX), dtype=np.float64)
    half = source.shape[1] // 2
    for j in range(NY):
        rows = np.asarray(source[j * 20 : (j + 1) * 20], dtype=np.float64)
        transformed = np.concatenate((rows[:, half:], rows[:, :half]), axis=1)[:, ::-1]
        result[j] = transformed.reshape(20, NX, 20).mean(axis=(0, 2))
    return result


def mean_mola64_patch(path: Path) -> np.ndarray:
    """Aggregate one 180-by-90 degree MOLA quadrant from 1/64 to 1 degree."""
    require_size(path, 11520 * 5760 * 2)
    source = np.memmap(path, dtype=">i2", mode="r", shape=(5760, 11520))
    result = np.empty((90, 180), dtype=np.float64)
    for j in range(90):
        y0 = 5760 - (j + 1) * 64
        y1 = 5760 - j * 64
        rows = np.asarray(source[y0:y1], dtype=np.float64)
        result[j] = rows.reshape(64, 180, 64).mean(axis=(0, 2))
    return result


def tes_to_wrf_longitudes(field: np.ndarray) -> np.ndarray:
    """Interpolate TES IAU-1994 bin centers onto 0.5..359.5 IAU-2000 centers."""
    source_lon = (TES_WRF_LON0 + np.arange(NX)) % 360.0
    order = np.argsort(source_lon)
    sorted_lon = source_lon[order]
    target_lon = np.arange(NX, dtype=np.float64) + 0.5
    result = np.empty_like(field, dtype=np.float64)
    for j, row in enumerate(field):
        values = row[order]
        extended_lon = np.concatenate(([sorted_lon[-1] - 360.0], sorted_lon, [sorted_lon[0] + 360.0]))
        extended_values = np.concatenate(([values[-1]], values, [values[0]]))
        result[j] = np.interp(target_lon, extended_lon, extended_values)
    return result


def add_border(field: np.ndarray) -> np.ndarray:
    """Add periodic longitude and replicated polar borders expected by geogrid."""
    if field.ndim == 2:
        field = field[np.newaxis, :, :]
    if field.shape[1:] != (NY, NX):
        raise ValueError(f"unexpected tile shape {field.shape}")
    with_y_border = np.pad(field, ((0, 0), (BORDER, BORDER), (0, 0)), mode="edge")
    return np.concatenate(
        (with_y_border[:, :, -BORDER:], with_y_border, with_y_border[:, :, :BORDER]),
        axis=2,
    )


def write_index(directory: Path, lines: list[str]) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    common = [
        "projection=regular_ll",
        "dx=1.0",
        "dy=1.0",
        "known_x=1.0",
        "known_y=1.0",
        "known_lat=-89.5",
        "known_lon=0.5",
        f"tile_x={NX}",
        f"tile_y={NY}",
        f"tile_bdr={BORDER}",
    ]
    (directory / "index").write_text("\n".join(lines + common) + "\n", encoding="ascii")


def write_tile(directory: Path, field: np.ndarray, dtype: str) -> dict[str, object]:
    encoded = add_border(field).astype(np.dtype(dtype), copy=False)
    path = directory / TILE_NAME
    encoded.tofile(path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return {
        "path": f"{directory.name}/{path.name}",
        "bytes": path.stat().st_size,
        "sha256": digest,
        "levels": int(encoded.shape[0]),
    }


def build(data_root: Path, output_root: Path) -> None:
    mola_dir = data_root / "topo/topo_latlon/topo64"
    albedo_new_path = data_root / "albedo/AlbMY24map.bin"
    albedo_old_path = data_root / "old_surface/tes_albedo_filled.dat"
    inertia_path = data_root / "old_surface/tes_inertia_filled.dat"
    roughness_path = data_root / "old_surface/roughness.dat"
    soil_ice_pc_path = data_root / "subsurface/grs_subsurface_ice.dat"
    soil_ice_dp_path = data_root / "subsurface/h2odepth.dat"

    require_size(albedo_new_path, 14400 + 7200 * 3600 * 2)
    check_fortran_record(albedo_old_path, 3040 * 1440 * 4)
    check_fortran_record(inertia_path, 7600 * 3600 * 4)
    check_fortran_record(roughness_path, 720 * 360)
    check_fortran_record(soil_ice_pc_path, 72 * 36 * 4)
    check_fortran_record(soil_ice_dp_path, 360 * 180 * 4)

    # topo64=.true. selects four 1/64-degree MOLA quadrants. Each source patch
    # is north-to-south; assemble a south-to-north, 0..360 east WPS grid.
    mola_ne = mean_mola64_patch(mola_dir / "megt90n000gb.img")
    mola_nw = mean_mola64_patch(mola_dir / "megt90n180gb.img")
    mola_se = mean_mola64_patch(mola_dir / "megt00n000gb.img")
    mola_sw = mean_mola64_patch(mola_dir / "megt00n180gb.img")
    topography = np.block([[mola_se, mola_sw], [mola_ne, mola_nw]])

    # MarsWRF first loads the old filled albedo globally, then overwrites |lat| < 87
    # with the MY24 product selected by alb_my=1.
    old_albedo_raw = np.memmap(
        albedo_old_path, dtype=">f4", mode="r", offset=4, shape=(1440, 3040)
    )
    old_albedo = mean_blocks(old_albedo_raw[:, 80:2960], 8, 8)
    my24_raw = np.memmap(
        albedo_new_path, dtype=">i2", mode="r", offset=14400, shape=(3600, 7200)
    )
    my24_albedo = mean_tes_map(my24_raw) / 10000.0
    albedo_natural = old_albedo.copy()
    lat_centers = np.arange(NY, dtype=np.float64) - 89.5
    use_my24 = (lat_centers > -87.0) & (lat_centers < 87.0)
    albedo_natural[use_my24] = my24_albedo[use_my24]
    albedo = tes_to_wrf_longitudes(albedo_natural)

    # ti2007=0 in em_global_mars/namelist.input, so use the old filled TES map.
    inertia_raw = np.memmap(
        inertia_path, dtype=">f4", mode="r", offset=4, shape=(3600, 7600)
    )
    inertia_core = np.maximum(inertia_raw[:, 200:7400], 10.0)
    thermal_inertia = tes_to_wrf_longitudes(mean_blocks(inertia_core, 20, 20))

    # MOLA pulse width is a signed byte in the file; MarsWRF maps it to 0..255
    # and then to roughness length with 0.15 / 198 m per count.
    roughness_raw = np.memmap(
        roughness_path, dtype="i1", mode="r", offset=4, shape=(360, 720)
    )
    pulse_width = np.asarray(roughness_raw, dtype=np.int16)
    pulse_width[pulse_width < 0] += 256
    roughness = mean_blocks(pulse_width, 2, 2) * (0.15 / 198.0)
    roughness = np.roll(roughness, -180, axis=1)

    # The ideal Mars setup reads these fields directly. The 5-degree GRS map
    # is expanded without inventing sub-grid variation; depth is already 1 degree.
    soil_ice_pc_raw = np.memmap(
        soil_ice_pc_path, dtype=">f4", mode="r", offset=4, shape=(36, 72)
    )
    soil_ice_pc = np.repeat(np.repeat(np.asarray(soil_ice_pc_raw), 5, axis=0), 5, axis=1)
    soil_ice_pc = np.roll(soil_ice_pc, -180, axis=1)
    soil_ice_dp_raw = np.memmap(
        soil_ice_dp_path, dtype=">f4", mode="r", offset=4, shape=(180, 360)
    )
    soil_ice_dp = np.roll(np.asarray(soil_ice_dp_raw), -180, axis=1)

    if not (-9000.0 < float(topography.min()) < float(topography.max()) < 22000.0):
        raise ValueError("MOLA values are outside the expected Mars elevation range")
    if not (0.0 < float(albedo.min()) < float(albedo.max()) < 1.0):
        raise ValueError("albedo values are outside 0..1")
    if not (10.0 <= float(thermal_inertia.min()) < float(thermal_inertia.max()) < 10000.0):
        raise ValueError("thermal inertia values are outside the expected range")

    output_root.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, object]] = []

    def emit(
        name: str,
        field: np.ndarray,
        dtype: str,
        index_lines: list[str],
        scale_factor: float = 1.0,
    ) -> None:
        directory = output_root / name
        write_index(directory, index_lines)
        record = write_tile(directory, field, dtype)
        record.update(
            {
                "dataset": name,
                "encoded_minimum": float(np.min(field)),
                "encoded_maximum": float(np.max(field)),
                "decoded_minimum": float(np.min(field)) * scale_factor,
                "decoded_maximum": float(np.max(field)) * scale_factor,
            }
        )
        records.append(record)

    albedo_i2 = np.rint(np.clip(albedo, 0.0, 1.0) * 10000.0).astype(np.uint16)
    inertia_i4 = np.rint(thermal_inertia).astype(np.int32)
    roughness_i2 = np.rint(np.clip(roughness, 0.0, None) * 100000.0).astype(np.uint16)

    emit(
        "mola_topography",
        np.rint(topography).astype(np.int16),
        ">i2",
        [
            "type=continuous",
            "signed=yes",
            "wordsize=2",
            "endian=big",
            "tile_z=1",
            "scale_factor=1.0",
            'units="m"',
            'description="1-degree block mean of topo64 MOLA selected by MarsWRF"',
        ],
    )
    emit(
        "mars_surface_class",
        np.ones((NY, NX), dtype=np.uint8),
        "u1",
        [
            "type=categorical",
            "signed=no",
            "wordsize=1",
            "tile_z=1",
            "category_min=1",
            "category_max=24",
            'units="category"',
            'description="Mars all-land class in a 24-category compatibility axis"',
        ],
    )
    emit(
        "mars_soil_class",
        np.ones((NY, NX), dtype=np.uint8),
        "u1",
        [
            "type=categorical",
            "signed=no",
            "wordsize=1",
            "tile_z=1",
            "category_min=1",
            "category_max=16",
            'units="category"',
            'description="Single neutral Mars soil class in a 16-category axis"',
        ],
    )
    emit(
        "mars_deep_soil_temp",
        np.full((NY, NX), 200, dtype=np.int16),
        ">i2",
        [
            "type=continuous",
            "signed=yes",
            "wordsize=2",
            "endian=big",
            "tile_z=1",
            "scale_factor=1.0",
            'units="K"',
            'description="200 K compatibility field; Data has no deep-soil climatology"',
        ],
    )
    emit(
        "mars_greenfrac",
        np.zeros((12, NY, NX), dtype=np.uint8),
        "u1",
        [
            "type=continuous",
            "signed=no",
            "wordsize=1",
            "tile_z_start=1",
            "tile_z_end=12",
            "scale_factor=0.01",
            'units="fraction"',
            'description="Twelve zero vegetation layers for Mars"',
        ],
        scale_factor=0.01,
    )
    emit(
        "mars_albedo12m",
        np.repeat(albedo_i2[np.newaxis, :, :], 12, axis=0),
        ">u2",
        [
            "type=continuous",
            "signed=no",
            "wordsize=2",
            "endian=big",
            "tile_z_start=1",
            "tile_z_end=12",
            "scale_factor=0.01",
            'units="percent"',
            'description="MY24 TES albedo repeated monthly; geogrid values are percent"',
        ],
        scale_factor=0.01,
    )
    emit(
        "mars_snow_albedo",
        np.zeros((NY, NX), dtype=np.uint8),
        "u1",
        [
            "type=continuous",
            "signed=no",
            "wordsize=1",
            "tile_z=1",
            "scale_factor=1.0",
            'units="percent"',
            'description="Zero compatibility field; Data has no snow-albedo climatology"',
        ],
    )
    emit(
        "mars_albedo",
        albedo_i2,
        ">u2",
        [
            "type=continuous",
            "signed=no",
            "wordsize=2",
            "endian=big",
            "tile_z=1",
            "scale_factor=0.0001",
            'units="fraction"',
            'description="MarsWRF alb_my=1 surface albedo reconstructed from Data"',
        ],
        scale_factor=0.0001,
    )
    emit(
        "mars_thermal_inertia",
        inertia_i4,
        ">i4",
        [
            "type=continuous",
            "signed=yes",
            "wordsize=4",
            "endian=big",
            "tile_z=1",
            "scale_factor=1.0",
            'units="J m-2 K-1 s-0.5"',
            'description="MarsWRF ti2007=0 TES thermal inertia reconstructed from Data"',
        ],
    )
    emit(
        "mars_emissivity",
        np.full((NY, NX), 10000, dtype=np.uint16),
        ">u2",
        [
            "type=continuous",
            "signed=no",
            "wordsize=2",
            "endian=big",
            "tile_z=1",
            "scale_factor=0.0001",
            'units="fraction"',
            'description="Emissivity 1.0 used by MarsWRF ideal surface setup"',
        ],
        scale_factor=0.0001,
    )
    emit(
        "mars_h2oice",
        np.zeros((NY, NX), dtype=np.int32),
        ">i4",
        [
            "type=continuous",
            "signed=yes",
            "wordsize=4",
            "endian=big",
            "tile_z=1",
            "scale_factor=1.0",
            'units="kg m-2"',
            'description="Zero initial H2O ice for the configured mp_physics=48"',
        ],
    )
    emit(
        "mars_roughness",
        roughness_i2,
        ">u2",
        [
            "type=continuous",
            "signed=no",
            "wordsize=2",
            "endian=big",
            "tile_z=1",
            "scale_factor=0.00001",
            'units="m"',
            'description="MOLA pulse-width roughness using MarsWRF 0.15/198 conversion"',
        ],
        scale_factor=0.00001,
    )
    emit(
        "mars_soil_ice_pc",
        np.rint(soil_ice_pc * 10000.0).astype(np.uint16),
        ">u2",
        [
            "type=continuous",
            "signed=no",
            "wordsize=2",
            "endian=big",
            "tile_z=1",
            "scale_factor=0.0001",
            'units="fraction"',
            'description="GRS subsurface H2O ice volume fraction from MarsWRF Data"',
        ],
        scale_factor=0.0001,
    )
    emit(
        "mars_soil_ice_dp",
        np.rint(soil_ice_dp * 10000.0).astype(np.int32),
        ">i4",
        [
            "type=continuous",
            "signed=yes",
            "wordsize=4",
            "endian=big",
            "tile_z=1",
            "scale_factor=0.0001",
            "missing_value=-9999.0",
            'units="m"',
            'description="Depth to subsurface H2O ice from MarsWRF Data"',
        ],
        scale_factor=0.0001,
    )

    manifest = {
        "source_root": str(data_root.resolve()),
        "grid": {"nx": NX, "ny": NY, "dx_degrees": 1.0, "dy_degrees": 1.0},
        "marswrf_selection": {"topo64": True, "alb_my": 1, "ti2007": 0},
        "datasets": records,
    }
    (output_root / "build_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="ascii"
    )


def main() -> None:
    project = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root",
        type=Path,
        default=project.parents[1] / "DATA/Data",
        help="MarsWRF Data directory",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=project / "sample/WPS_GEOG_MARS",
        help="output WPS geogrid data directory",
    )
    args = parser.parse_args()
    build(args.data_root, args.output_root)


if __name__ == "__main__":
    main()
