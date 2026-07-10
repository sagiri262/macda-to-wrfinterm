# MACDA v2.0 to MarsWRF intermediate

The MACDA path is independent of the retained Earth/CMIP converter. Install its
small dependency set and run it from this directory:

```bash
pip install -r requirements.macda.txt
python run_macda2w.py --dry-run
python run_macda2w.py --time-index 0
python run_macda2w.py --max-times 2
```

The default input is the MY28/SOY507-SOY537 file under
`/home/zy/WRF/DATA/macda-mro-mcs`. Output goes to `output/`.

## Time conversion

MACDA has three independent time descriptions: continuous `time` in sols,
`MY_Ls/Ls`, and a `Mars_date` sol calendar. MarsWRF's checked files use
`YYYY-DDDDD_HH:MM:SS`, where `DDDDD` is a 1-based sol of year. The converter
uses the documented MACDA five-year cycle (`669,668,669,668,669`) and month
lengths to make that value, then checks it against continuous `time`.

```text
+0028-10-07T02:00:00A -> 0028-00507_02:00:00.0000
```

The old direct replacement produced `0028-10-07_...`, which is not a valid
MarsWRF planetary date. `XFCST` remains elapsed Martian hours from the first
selected record.

## Field pipeline

The active field table is `db/MACDA-v2_CORE.csv`:

```text
temp    -> TT          sigma -> fixed pressure
uwind   -> UU          sigma -> fixed pressure, geographic-relative wind
vwind   -> VV          sigma -> fixed pressure, geographic-relative wind
derived -> PRESSURE    metgrid derives final PRES
geop    -> GHT         geopotential / Mars gravity
psurf   -> PSFC
tsurf   -> SKINTEMP
zero    -> SPECHUMD, QV
coldust -> TAU_OD2D    normalized to 700 Pa, optional
co2ice  -> CO2ICE      optional
```

Pressure is `p(k,j,i) = psurf(j,i) * lev(k)` and vertical interpolation is
linear in log-pressure. Each written file is immediately read back to verify
HDATE, required field names, dimensions, finite values, byte order, and wind
flags.

`omega`, `swflux`, and `lwflux` are diagnostics rather than real.exe initial
state. `dustmmr` cannot be assigned to MarsWRF's two dust bins without a
scientific bin-partition assumption, so these variables are intentionally not
mislabelled as WRF state.

Optional `TAU_OD2D` and `CO2ICE` require the entries in
`sample/MACDA-v2/METGRID.TBL.MARS.additions` to be added to the active Mars
METGRID table. See `MARSWRF_AUDIT.md` for the source-level time and physics
audit and for downstream limitations.
