# 0004 HDATE labels across a Mars year boundary

## Added

- `macda2wrf/mars_time.py`: `WrfPlanetDate`, `MARSWRF_FIXED_YEAR_SOLS`,
  `HDATE_STRATEGIES`, `whole_sols_since_macda_epoch`, `relabel_fixed_year`.
- `hdate_strategy = marswrf_fixed669`, which relabels records onto the fixed
  669-sol calendar, anchored at the first record of the run.
- `FixedYearLabelTest` in `tests/test_mars_time.py`: anchor, boundary,
  post-boundary offset, single-year equivalence, unknown strategy name.
- The cause and fix sections of `docs/跨年模拟`.

## Changed

- `macda2wrf/converter.py`: `_write_one_time` takes the `Mars_date` of the
  first record of the run as the anchor, prints one line whenever a label
  differs from the true MACDA date, and `--dry-run` now reports
  `hdate_strategy` plus the first and last labels.
- `macda2wrf/config.py`: validates `hdate_strategy` against
  `HDATE_STRATEGIES` at load time instead of failing while writing.
- `macda2wrf/mars_time.py`: `MarsDate.wrf_date`/`hdate`/`filename_stamp`
  delegate to `WrfPlanetDate`; `sols_since_macda_epoch` is now the whole-sol
  count plus the fraction of sol.
- `config/config.MACDA-v2-MY34-MY35.ini`: uses `marswrf_fixed669` and clears
  `max_times` for a full conversion.
- `docs/MARSWRF_AUDIT.md`, `docs/MARSWRF_AUDIT.zh-CN.md`: record this
  no-Fortran-change option after the cross-year limitation.

## Removed

No files were removed.

## Design rationale

metgrid steps dates on a fixed 669-sol Mars year
(`WPS/metgrid/src/module_date_pack.F:17`), while MY34 has 668 sols in the MACDA
sol calendar. Across the boundary the converter wrote `MACDA:0035-00001_00`
but metgrid asked for `MACDA:0034-00669_00`, which produced
`Couldn't open file` followed by `mandatory field TT was not found`.

Patching the WPS and WRF calendars together is the most faithful fix, but 669
sols is also hard-coded in places such as `WRF/external/esmf_time_f90/Meat.F90`,
so patching WPS alone would make `real.exe`/`wrf.exe` fail at the same boundary
in the same way — the largest effort and the largest regression risk. Matching
MarsWRF's fixed year length on the converter side instead keeps sol spacing
identical to `geth_newdate`/`geth_idts`, so every file metgrid asks for exists,
and no Fortran changes.

The price is that labels after the boundary lag the true MACDA date by one sol:
the extra sol WPS believes in becomes a phantom MY34 sol 669. The time axis
stays strictly uniform, but comparisons against observations must subtract that
sol, and the model's internal L_s runs about 0.5 degrees early. The strategy is
therefore an explicit config option; the default stays `marswrf_sol` and no
single-year run changes behaviour.

## Verification

```text
$ python -m unittest discover -s tests -t .
Ran 14 tests in 0.021s
OK
```

The four time steps around the boundary of the MY34SOY665_MY35SOY027 file:

```text
[macda2w] marswrf_fixed669 relabel 0035-00001_00 -> 0034-00669_00
[macda2w] marswrf_fixed669 relabel 0035-00001_22 -> 0034-00669_22
[macda2w] marswrf_fixed669 relabel 0035-00002_00 -> 0035-00001_00
output/MY34-MY35/MACDA:0034-00668_22
output/MY34-MY35/MACDA:0034-00669_00
output/MY34-MY35/MACDA:0034-00669_22
output/MY34-MY35/MACDA:0035-00001_00
```

`--dry-run` reports 360 time steps with first and last labels
`0034-00665_02:00:00` and `0035-00026_00:00:00`, exactly 359 intervals of
7200 seconds on the fixed 669-sol calendar.

## Open issues

- `namelist.wps` still needs `end_date` changed to `0035-00026_00:00:00` by
  hand; the converter does not write namelists.
- The full 360-step conversion and `metgrid.exe` have not been rerun. Old
  `MACDA:*` and `met_em.*` files must be deleted from both the local output
  directory and the cluster WPS directory first, or files that do not belong
  to the new sequence will be mixed in.
- The one-sol label offset after the boundary has to be carried into the
  InSight comparison. Removing it entirely still requires patching the WPS and
  WRF calendars to the 669/668 cycle as described in section 5 of
  `docs/MARSWRF_AUDIT.md`, or a correctly epoch-aligned MARS24/MSD build.
