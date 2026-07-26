"""Synthetic MACDA records around the MY 34 to MY 35 year boundary.

The real forcing file, ``mro-mcs-reanalysis_mars_MY34SOY665_MY35SOY027_v2-0.nc``,
is 806 MB and lives outside this repository. These helpers rebuild its time axis
from the MACDA calendar itself, so the unit tests never read that file, the
``output/`` tree, or ``public/``. ``tests.test_mars_time`` asserts that the
rebuilt endpoints match the real file.
"""

from __future__ import annotations

from typing import Sequence

from macda2wrf.mars_time import MacdaRecord, mars_date_at, whole_sols_since_epoch


FIRST_YEAR = 34
FIRST_SOL_OF_YEAR = 665
FIRST_HOUR = 2
STEP_MARTIAN_HOURS = 2
RECORD_COUNT = 360

# Selected time indices of the real file: the last true MY 34 record, the first
# three true MY 35 records, and the final record.
LAST_MY34_INDEX = 46
FIRST_MY35_INDEX = 47
LAST_INDEX = RECORD_COUNT - 1


def boundary_records(count: int = RECORD_COUNT) -> list[MacdaRecord]:
    """Return ``count`` records two Martian hours apart from MY 34 sol 665 02:00."""

    if count < 1:
        raise ValueError(f"count must be at least 1, got {count}")
    epoch_sols = whole_sols_since_epoch(FIRST_YEAR, FIRST_SOL_OF_YEAR)
    records: list[MacdaRecord] = []
    for offset in range(count):
        sol_offset, hour = divmod(FIRST_HOUR + STEP_MARTIAN_HOURS * offset, 24)
        mars_date = mars_date_at(epoch_sols + sol_offset, hour=hour)
        records.append(MacdaRecord(time_index=offset, mars_date=mars_date.macda_string,
                                   time_sols=epoch_sols + sol_offset + hour / 24.0))
    return records


def time_axis(records: Sequence[MacdaRecord]) -> list[float]:
    """Return the MACDA ``time`` coordinate indexed by time index."""

    axis = [0.0] * (max(record.time_index for record in records) + 1)
    for record in records:
        axis[record.time_index] = record.time_sols
    return axis


class FakeMacdaReader:
    """Stand-in for :class:`macda2wrf.macda_reader.MacdaReader`.

    It serves only the ``Mars_date`` strings the labelling path reads, so no
    NetCDF4/HDF5 file and no h5py/netCDF4 install is needed.
    """

    def __init__(self, records: Sequence[MacdaRecord]):
        self.mars_dates = {record.time_index: record.mars_date for record in records}
        self.backend = "fixture"

    def read_time_string(self, time_index: int) -> str:
        if time_index not in self.mars_dates:
            raise KeyError(f"fixture has no record at time index {time_index}")
        return self.mars_dates[time_index]
