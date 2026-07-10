# MarsWRF time, physics, and MACDA conversion audit

This audit is based on the checked source tree and the actual files below:

- `marswrf/WRFV3/share` and `marswrf/WRFV3/external/esmf_time_f90`
- `marswrf/WRFV3/phys`
- `DATA/macda-mro-mcs/mro-mcs-reanalysis_mars_MY28SOY507_MY28SOY537_v2-0.nc`
- `marswrf/em_global_mars/wrfinput_d01`
- `marswrf/em_global_mars/wrfout_d01_0003-00001_00:00:00`
- `marswrf/em_global_mars/wrfrst/wrfrst_d01_0001-00520_00:00:00`

## 1. How MarsWRF represents time

MarsWRF has two compile-time time modes. `configure` selects `mars` for the
traditional mode and `mars24`/`marsda` for `MARS24_TIMING`.

### Traditional mode used by the checked files

`WRF_PLANET` changes the WRF string layout from Gregorian
`YYYY-MM-DD_HH:MM:SS` to planetary `YYYY-DDDDD_HH:MM:SS`:

- `share/module_date_time.F:203-210, 228-235` parses columns 6-10 as one
  five-digit day/sol field and sets month to zero.
- `share/module_date_time.F:251-264, 541-560` validates the sol against
  `PLANET_YEAR`.
- `share/module_date_time.F:668-722, 742-765` rolls the sol across a planetary
  year and formats it with `I5.5`.
- `external/esmf_time_f90/Meat.F90:202-209` sets a Mars year to 669 model days.
- `external/esmf_time_f90/ESMF_Time.F90:751-756` emits the same five-digit
  planetary string.

One model day is represented by 86400 clock units so existing WRF alarms and
date arithmetic remain usable. This does not assert that a physical sol lasts
86400 SI seconds. `frame/module_driver_constants.F:95-101` defines
`P2SI=1.027491252`, or about 88775.2 SI seconds per sol, and the radiation and
surface physics apply that scale where SI time is required.

`share/set_timekeeping.F` constructs WRFU clocks from namelist year/day/hour
values and attaches history, restart, boundary, and auxiliary alarms.
`module_bc_time_utilities.F` only compares those clock objects to decide when
to read lateral boundaries. `wrf_timeseries.F` writes station diagnostics and
does not perform Earth/Mars conversion.

### MARS24 mode

With `MARS24_TIMING`, `PLANET_YEAR` is enlarged to 100000 so the five-digit day
field can carry a continuous Mars Solar Date (MSD), rather than a repeating
669-sol year. `share/module_mars24.F` implements:

- UTC Julian day to TT, including the leap-second offset (`:25-70`)
- TT/J2000 offset (`:72-77`)
- J2000 to areocentric solar longitude Ls (`:79-150`)
- MSD and J2000 conversion (`:167-179`)
- MTC/local mean/local true solar time (`:181-225`)

`share/module_planet_utilities.F:4-19` and
`phys/module_radiation_driver.F:2513-2630` use these routines for Ls,
declination, heliocentric distance, and equation of time. This is the only
source path that converts an Earth UTC/J2000 representation astronomically.
The normal NetCDF I/O path itself does not store an Earth date alongside the
Mars date.

### Actual WRF file evidence

The checked files are traditional mode, not MARS24 mode:

| File | `Times` | Key attributes |
|---|---|---|
| `wrfinput_d01` | `0001-00001_00:00:00` | `PLANET_YEAR=669`, `P2SI=1.027491` |
| `wrfout_d01_0003-00001_00:00:00` | `0003-00001_00:00:00` | same |
| `wrfrst_d01_0001-00520_00:00:00` | `0001-00520_00:00:00` | same |

`share/output_wrf.F:709-762` writes the planetary constants to NetCDF.
`share/input_wrf.F:174-217` reads `SIMULATION_START_DATE` with the planetary
year/sol layout. `START_DATE`, `SIMULATION_START_DATE`, the `Times` variable,
and filenames therefore all use the same MarsWRF string convention.

## 2. MACDA time alignment

The input file has 360 records at two Martian-hour cadence. Its independent
coordinates state:

- `time=3180.08333333333 .. 3210` sols since MY24 sol 1 at 00:00 MTC
- `MY_Ls=28`
- `Ls=264.283 .. 283.5132` degrees
- `Mars_date=+0028-10-07T02:00:00A .. +0028-10-37T00:00:00A`

The documented MACDA sol calendar repeats year lengths
`669,668,669,668,669`. Its month lengths in a 668-sol year are:

```text
56, 55, 56, 55, 56, 56, 55, 56, 55, 56, 56, 56
```

Month 12 has 57 sols in a 669-sol year. The first nine months contain 500
sols, so MY28 month 10 sol 7 is sol-of-year 507. The conversion is therefore:

```text
+0028-10-07T02:00:00A
  -> MY 28, SOY 507, 02:00:00 MTC
  -> MarsWRF 0028-00507_02:00:00
  -> HDATE   0028-00507_02:00:00.0000
```

The converter independently reconstructs continuous time:

```text
MY24 668 + MY25 669 + MY26 669 + MY27 668
+ (SOY507 - 1) + 2/24 = 3180.08333333333 sols
```

A mismatch between this value and the NetCDF `time` coordinate is now fatal.
`XFCST` uses the difference in continuous `time` multiplied by 24 Martian
hours, so it remains continuous within the selected forcing sequence.

The old converter merely removed `+`, `A`, and `T`; it would have passed
`0028-10-07_02:00:00` to software that expects columns 6-10 to be one sol
number. That was not a MarsWRF date.

## 3. MACDA to WPS/WRF fields

The active mapping and transformations are:

| MACDA | Intermediate | Processing | Downstream role |
|---|---|---|---|
| `temp` | `TT` | sigma to pressure | final `TT` |
| `uwind` | `UU` | sigma to pressure | final staggered `UU` |
| `vwind` | `VV` | sigma to pressure | final staggered `VV` |
| derived | `PRESSURE` | constant slab per pressure level | metgrid derives `PRES` |
| `geop` | `GHT` | divide by 3.72 m/s2 | final `GHT` |
| `psurf` | `PSFC` | identity, Pa | surface pressure |
| `tsurf` | `SKINTEMP` | identity, K | skin temperature |
| zero | `SPECHUMD` | pressure-level zero | dry initial state |
| zero | `QV` | pressure-level zero | dry water-vapor state |
| `coldust` | `TAU_OD2D` | `coldust/psurf*700 Pa` | optional Mars dust opacity |
| `co2ice` | `CO2ICE` | identity | optional surface CO2 ice |

MACDA pressure is exactly `p(k,j,i)=psurf(j,i)*lev(k)`. The 3-D fields are
interpolated linearly in log-pressure to configured Mars levels. Latitude is
reversed from north-to-south source order, longitude is normalized from
`[-180,175]` to `[0,355]`, and winds are marked geographic-relative.

`omega` is pressure vertical velocity and cannot be copied into WRF geometric
`W`. `swflux` and `lwflux` are diagnostic surface fluxes that MarsWRF
recomputes. `dustmmr` is one 1.5-micron distribution, while this MarsWRF setup
uses two dust tracers/bins; assigning it to `TRC01/TRC02` requires an explicit
scientific partition and is intentionally not guessed.

## 4. Mars-related physics files

Dedicated Mars implementations:

- `module_ra_mars_common.F`: dust profiles, MCD/MGS/Viking/MCS/TES dust,
  aerosol shortwave and longwave heating (`:47-1449`)
- `module_ra_mars_kdm.F`: correlated-k Mars radiation (`:142-2944`)
- `module_ra_mars_wbm.F`: WBM visible/IR radiation (`:15-1005`)
- `module_ra_mars_burk.F`: Burke shortwave heating (`:11-377`)
- `module_ra_mars_uv.F`: UV heating (`:110-494`)
- `module_ra_houben.F`: simplified Mars Newtonian cooling (`:28-145`)
- `module_mp_mars_common.F`: Mars particle sedimentation/Stokes support
  (`:18-454`)
- `module_mp_mars_co2_simple.F`: CO2 condensation/sublimation (`:11-305`)
- `module_mp_mars_h2o_simple.F`: H2O cloud microphysics (`:11-164`)
- `module_mp_mars_basudustlifting.F`: Basu dust lifting (`:13-234`)
- `module_mp_sedim_dust.F`, `module_mp_sedim_water.F`: dust/water settling
- `module_sf_mars_cendustlifting.F`: two-bin, one-bin, variable-tau, and
  prescribed dust injection (`:22-787`)
- `module_mp_mars_chem.F`: passive chemical tracer tendencies (`:14-408`)
- `module_sf_planet_simple.F`: Mars surface/subsurface energy solver and
  initialization (`:12-812`, especially the Mars selection at `:416-530`)

Shared drivers with substantial `WRF_MARS` blocks:

- `module_radiation_driver.F`: selects Mars radiation, orbital/Ls calculation,
  dust and cloud optical profiles (`:119-209, 616-670, 1434-1541,
  1875-2052, 2334-2408, 2513-2946`)
- `module_microphysics_driver.F`: connects CO2, H2O, dust sedimentation,
  lifting, and tracers (`:52-142, 309-386, 466-509, 1149-1543`)
- `module_surface_driver.F`: passes Mars ice/dust/thermal-inertia fields and
  invokes planetary ground physics (`:126-212, 530-560, 1160-1301,
  2148-2220, 3662-3810`)
- `module_pbl_driver.F`: passes Mars tracers and planetary boundary-layer
  arguments (`:77-136, 570-688, 871-878, 1067-1148, 1244-1251, 1591-1598`)
- `module_physics_init.F`: initializes Mars radiation, surface, dust, CO2/H2O,
  and tracer schemes (`:33-38, 253-259, 572-633, 1085-1093,
  1289-1502, 1647-1654, 1997-2077`)
- `module_physics_addtendc.F`: couples Mars scalar/tracer tendencies into the
  dynamics (`:40-44, 93-100, 145-149, 237-282, 365-414, 479-528, 695-744`)
- `module_bl_mrf.F`, `module_bl_myjpbl.F`, `module_bl_ysu.F`: CO2/Mars
  thermodynamics and extra tracer transport
- `module_sf_sfclay.F` and `_3012.F`: Mars gas constants and surface-layer
  treatment

## 5. Verified scope and remaining downstream gap

The Python path now completes and validates MACDA to WRF intermediate output.
WPS 4.6 `rd_intermediate.exe` also reads all 137 records successfully and
reports the expected 72 x 36 Mars grid, 3389.92 km radius, source `MACDA`, and
planetary HDATE.
The checked `marswrf/WPS` date utilities do not contain the `WRF_PLANET`
five-digit-sol branches present in WRFV3, and there is no local `metgrid.exe`
or `real.exe` to execute an end-to-end WPS/real test. Consequently the correct
Mars HDATE may still require porting the WRF planetary date logic into this old
WPS tree before `metgrid` can schedule it. The sample namelist records the
correct target dates but does not conceal that source-level limitation.

The optional `TAU_OD2D` and `CO2ICE` records also need the supplied
`METGRID.TBL.MARS.additions` entries. A Mars static geography dataset remains
required for geogrid/real.

Finally, traditional MarsWRF always treats every year as 669 sols, while the
MACDA sol calendar includes 668-sol years. The requested MY28 file does not
cross a year boundary. A forcing run that crosses one should either patch the
WRF/WPS calendar consistently or use a correctly epoch-aligned MARS24/MSD
build; it should not silently rely on the traditional year arithmetic.
