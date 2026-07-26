# 0005 Stable label anchor for cross-Mars-year runs

## Added

- `macda2wrf/mars_time.py`: `MacdaRecord`, `LabelledRecord`, `make_label`,
  `build_label_sequence`, `check_unique_labels`, `check_label_spacing`,
  `check_hdate_strategy`, `label_sols`, `fixed_year_sols`, `mars_date_at`,
  `month_and_sol`, `whole_sols_since_epoch`, `short_mars_years_left`,
  `ANCHORED_HDATE_STRATEGIES`, `MarsDate.macda_string`,
  `WrfPlanetDate.sol_fraction`.
- `macda2wrf/converter.py`: `_anchor_index`, `_select_indices`,
  `_validate_indices`, `_build_labels`.
- `tests/macda_boundary_fixture.py`: rebuilds the time axis of
  `mro-mcs-reanalysis_mars_MY34SOY665_MY35SOY027_v2-0.nc` from the MACDA
  calendar, so no test reads that 806 MB file, `output/`, or `public/`.
- `tests/test_mars_time.py`: `MacdaCalendarInverseTest`,
  `StableAnchorLabelTest`, `LabelSequenceValidationTest`, and two more cases in
  `FixedYearLabelTest` (missing anchor, default strategy keeps true MY35).
- `tests/test_converter.py`: `ConverterStableAnchorTest` and four more cases in
  `ConverterIndexTest`.
- `docs/MACDA跨年固定669兼容方案.md`.

## Changed

- `macda2wrf/mars_time.py`: `marswrf_fixed669` now requires an explicit anchor.
  `make_hdate` no longer treats the record as its own anchor when the anchor is
  missing; it raises instead. `relabel_fixed_year` documents that it is a
  compatibility relabelling only.
- `macda2wrf/converter.py`: the label and `XFCST` anchor is the record at the
  configured `start_index`, read once per run, instead of the first
  `time_index` of the current command. Explicit `--time-index` values are
  sorted and deduplicated, rejected when outside the file, and rejected when
  earlier than the anchor. The whole selection is labelled and validated in one
  place before any file is written. `--dry-run` validates and prints every
  selected record instead of only the first and last. A selection that leaves a
  MACDA year shorter than 669 sols while using true MACDA labels now warns.
- `config/config.MACDA-v2-MY34-MY35.ini`: comments only. `start_index = 0` is
  documented as the run's stable anchor that must not change between a full run
  and a partial rerun, and `hdate_strategy` points at the new document.
- `docs/MARSWRF_AUDIT.md`, `docs/MARSWRF_AUDIT.zh-CN.md`: one sentence each,
  correcting "anchored at the first record of the run" to the configured
  `start_index` and pointing at the new document.

## Deleted

No files were removed.

## Reasoning

Record 0004 anchored the fixed 669-sol calendar at the first record of the
current run. That made the label of one NC record depend on which time indices
a single command happened to select:

```text
--time-index 47                  ->  0035-00001_00
--time-index 0 --time-index 47   ->  0034-00669_00
```

So a partial rerun, a resumed run, and a full run wrote different filenames for
the same data, and `XFCST` moved with them. The anchor has to be a property of
the dataset and the configuration, not of the command line, so it is now
`start_index`. A single record carries no information about where the fixed
calendar starts, which is why the anchor is mandatory rather than defaulted:
silently anchoring on the record itself is the bug, not a fallback.

Validation moved from per-file to per-sequence for the same reason. metgrid
derives filenames arithmetically, so what matters is not that each label is
well formed but that the sequence is unique (no output file overwrites
another), strictly increasing, and spaced exactly like the MACDA `time`
coordinate. `build_label_sequence` checks all three before the first file is
written, and `--dry-run` reports the whole selection.

## Validation

The tests were written but **not executed**: this change was made under a
local-only constraint, so no conversion, no metgrid, and no test run happened
here. `tests/test_mars_time.py` and `tests/test_converter.py` cover the full
run versus single-time-step equivalence, the MY34 sol 668 to fixed sol 669
mapping, label continuity across the boundary, `XFCST` relative to the config
anchor in a partial rerun, the missing-anchor error, and the unchanged default
`marswrf_sol` labels.

Static checks that were run locally: `ast.parse` on all five touched Python
files, no line over 120 characters, no tabs, and no bracket held open with a
single value per line.

The HPC verification checklist is section 7 of
`docs/MACDA跨年固定669兼容方案.md`.

## Remaining issues

- Nothing has been verified on the cluster. Until section 7 of the new document
  is complete, this change is "implemented, pending verification".
- `namelist.wps` still has `end_date = '0035-00027_00:00:00'`, which is the true
  MACDA date of the last record, not its fixed-669 label. It must become
  `0035-00026_00:00:00`: 718 h / 2 h = 359 intervals, so 360 time steps.
- `output/MY34-MY35` holds output written under the old, unstable anchor. It is
  a mixture of two labelling conventions and cannot be used as verification.
- The one-sol label offset after the boundary remains. Removing it still
  requires patching the WPS and WRF calendars to the 669/668 cycle, as
  described in section 5 of `docs/MARSWRF_AUDIT.md`.
- Identifiers stay `snake_case` to match the rest of the package, which departs
  from the lowerCamelCase rule in the shared coding standard. Renaming would be
  a repository-wide refactor outside this change.
