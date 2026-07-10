# 0002 MarsWRF time and field alignment

## Changed

- Added the documented MACDA 669/668 sol calendar conversion and continuous
  `time` cross-check.
- Changed intermediate pressure input from `PRES` to `PRESSURE`; metgrid
  derives final `PRES` from it.
- Marked `UU` and `VV` as geographic-relative winds.
- Added NetCDF packed/missing-value handling and stricter configuration/table
  validation.
- Added WRF intermediate read-back validation and unit tests.
- Replaced empty sample CSV placeholders and invalid sample namelist content.
- Removed the unused `db/MACDA_LEV.csv` placeholder with its duplicate header;
  `db/MACDA-v2_CORE.csv` is the sole active MACDA field table.
- Documented the MarsWRF time implementation, actual WRF/MACDA metadata,
  Mars-related physics blocks, and downstream WPS limitations.

## Checked conversion

```text
+0028-10-07T02:00:00A
time=3180.08333333333 sol
-> 0028-00507_02:00:00.0000
```
