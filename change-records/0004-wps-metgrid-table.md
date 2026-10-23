# 0004 WPS metgrid table

## Changed files

- `../WPS/metgrid/METGRID.TBL.ARW.copy`
- `../WPS/metgrid/METGRID.TBL`
- `sample/metgrid/METGRID.TBL.ARW.copy`
- `sample/metgrid/METGRID.TBL`
- `sample/namelist.wps.copy`

## Changes

- Preserved the original `METGRID.TBL.ARW`.
- Created the exact runtime filename required by `metgrid.exe`.
- Retained the existing core rules for `PRESSURE -> PRES`, `TT`, `UU`, `VV`,
  `GHT`, `PSFC`, `SKINTEMP`, `SPECHUMD`, and `QV`.
- Added interpolation entries for optional MACDA fields `TAU_OD2D` and
  `CO2ICE`.
- Set `opt_metgrid_tbl_path='./metgrid'` in the modified sample namelist.

No WPS source was changed by this table update, so recompilation is not
required for these changes alone.
