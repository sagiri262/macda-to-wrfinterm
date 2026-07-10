# 0001 Initial MACDA Converter Scaffold

## Added

- `.gitignore`
- `output/.gitkeep`
- `run_macda2w.py`
- `macda2wrf/__init__.py`
- `macda2wrf/config.py`
- `macda2wrf/grid.py`
- `macda2wrf/macda_reader.py`
- `macda2wrf/vertical.py`
- `macda2wrf/wrf_intermediate.py`
- `macda2wrf/converter.py`
- `conf/config.MACDA-v2.ini`
- `db/MACDA-v2_CORE.csv`
- `requirements.macda.txt`
- `README_MACDA.md`
- `change-records/README.md`
- `change-records/0001-initial-macda-scaffold.md`

## Modified

- No existing source file was modified in this step.

## Deleted

- No file was deleted in this step.

## Reasoning

The copied CMIP6 project remains intact. The MACDA work is isolated in a new
entry point and package because MACDA uses one NetCDF4/HDF5 file with sigma
levels, Mars dates, and Mars pressure magnitudes rather than CMIP6-style
one-variable files and Earth pressure levels.

The first field table targets the minimal MarsWRF real-data path:

```text
TT, UU, VV, PRES, GHT, PSFC, SKINTEMP, SPECHUMD, QV
```

It also writes optional Mars fields:

```text
TAU_OD2D, CO2ICE
```

Those optional fields may require downstream MarsWRF/WPS changes before they
are consumed by `real.exe` or `wrf.exe`.

## Validation

- `python3 -m py_compile` passed for the new entry point and package modules.
- A temporary writer test created a WRF intermediate record in `/tmp`.
- A small synthetic sigma-to-pressure interpolation test returned finite
  output.
- Config and CSV parsing found the MACDA file path, the variable table, 11
  field definitions, a 36 x 72 target grid, and 19 pressure levels.
- `python3 run_macda2w.py --dry-run` reaches the real MACDA open step, then
  stops because this Python environment lacks `h5py` and `netCDF4`.

## Remaining issues

- The local Python environment must provide `h5py` or `netCDF4`; the current
  default environment did not have either package at scaffold time.
- Mars-specific fields are not guaranteed to pass through current
  `METGRID.TBL` and `realonly` Registry logic.
- The `hdate_strategy=mars_date` path emits MACDA Mars calendar strings in WRF
  intermediate records. This should be checked against the MarsWRF WPS date
  parser before running a long production conversion.
