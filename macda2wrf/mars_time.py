"""MACDA sol-calendar conversion to the MarsWRF planetary date format."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence
import re


MACDA_EPOCH_YEAR = 24
MARS_SECONDS_PER_SOL = 24 * 60 * 60

# Montabone et al. sol calendar. Month 12 has one extra sol in a 669-sol year.
MONTH_LENGTHS_668 = (56, 55, 56, 55, 56, 56, 55, 56, 55, 56, 56, 56)

# Traditional MarsWRF and the patched WPS date pack use a fixed 669-sol year;
# see WPS metgrid/src/module_date_pack.F.
MARSWRF_FIXED_YEAR_SOLS = 669

FIXED_YEAR_STRATEGY = "marswrf_fixed669"

# Strategies that label records on the fixed 669-sol MarsWRF calendar rather than
# on the MACDA calendar. A single record carries no information about where the
# fixed calendar starts, so these strategies always need the anchor record.
ANCHORED_HDATE_STRATEGIES = frozenset({FIXED_YEAR_STRATEGY})
HDATE_STRATEGIES = frozenset({"marswrf_sol", "mars_date", FIXED_YEAR_STRATEGY})

_MARS_DATE_RE = re.compile(r"^(?P<sign>[+-])(?P<year>\d{4})-(?P<month>\d{2})-(?P<sol>\d{2})"
                           r"T(?P<hour>\d{2}):(?P<minute>\d{2}):(?P<second>\d{2})A$")


def check_hdate_strategy(strategy: str) -> str:
    """Reject an unknown ``hdate_strategy`` before any label is produced."""

    if strategy not in HDATE_STRATEGIES:
        raise ValueError(f"Unsupported hdate_strategy {strategy!r}; use one of "
                         f"{sorted(HDATE_STRATEGIES)}")
    return strategy


def mars_year_length(year: int) -> int:
    """Return the MACDA sol-calendar year length.

    The calendar repeats 669, 668, 669, 668, 669 sols in five-year cycles.
    The cycles containing this dataset start at MY 21, MY 26, and so on.
    """

    return 668 if year % 5 in {2, 4} else 669


def mars_month_lengths(year: int) -> tuple[int, ...]:
    lengths = list(MONTH_LENGTHS_668)
    if mars_year_length(year) == 669:
        lengths[-1] += 1
    return tuple(lengths)


def month_and_sol(year: int, sol_of_year: int) -> tuple[int, int]:
    """Split a MACDA sol of year into its month and its sol within that month."""

    if sol_of_year < 1:
        raise ValueError(f"Sol of year must be at least 1, got {sol_of_year}")
    remaining = sol_of_year
    for month, length in enumerate(mars_month_lengths(year), start=1):
        if remaining <= length:
            return month, remaining
        remaining -= length
    raise ValueError(f"Sol of year {sol_of_year} is outside MY {year}, which has "
                     f"{mars_year_length(year)} sols")


@dataclass(frozen=True)
class WrfPlanetDate:
    """A MarsWRF planetary date: four-digit year plus five-digit sol of year."""

    year: int
    sol_of_year: int
    hour: int
    minute: int
    second: int

    @property
    def wrf_date(self) -> str:
        """The 19-character MarsWRF date: YYYY-DDDDD_HH:MM:SS."""

        if not 0 <= self.year <= 9999:
            raise ValueError(f"MarsWRF four-digit year cannot represent MY {self.year}")
        if not 1 <= self.sol_of_year <= 99999:
            raise ValueError(f"MarsWRF five-digit sol field cannot represent sol "
                             f"{self.sol_of_year}")
        return (f"{self.year:04d}-{self.sol_of_year:05d}_"
                f"{self.hour:02d}:{self.minute:02d}:{self.second:02d}")

    @property
    def hdate(self) -> str:
        """The 24-character HDATE stored in WRF intermediate records."""

        return f"{self.wrf_date}.0000"

    @property
    def filename_stamp(self) -> str:
        return self.wrf_date[:13]

    @property
    def sol_fraction(self) -> float:
        seconds = self.hour * 3600 + self.minute * 60 + self.second
        return seconds / MARS_SECONDS_PER_SOL


@dataclass(frozen=True)
class MarsDate:
    year: int
    month: int
    sol: int
    hour: int
    minute: int
    second: int

    @property
    def sol_of_year(self) -> int:
        return sum(mars_month_lengths(self.year)[: self.month - 1]) + self.sol

    @property
    def sol_fraction(self) -> float:
        seconds = self.hour * 3600 + self.minute * 60 + self.second
        return seconds / MARS_SECONDS_PER_SOL

    @property
    def wrf_planet_date(self) -> WrfPlanetDate:
        return WrfPlanetDate(year=self.year, sol_of_year=self.sol_of_year,
                             hour=self.hour, minute=self.minute, second=self.second)

    @property
    def wrf_date(self) -> str:
        """The 19-character MarsWRF date: YYYY-DDDDD_HH:MM:SS."""

        return self.wrf_planet_date.wrf_date

    @property
    def hdate(self) -> str:
        """The 24-character HDATE stored in WRF intermediate records."""

        return self.wrf_planet_date.hdate

    @property
    def filename_stamp(self) -> str:
        return self.wrf_planet_date.filename_stamp

    @property
    def macda_string(self) -> str:
        """Render this date back into the MACDA ``Mars_date`` form."""

        return (f"+{self.year:04d}-{self.month:02d}-{self.sol:02d}"
                f"T{self.hour:02d}:{self.minute:02d}:{self.second:02d}A")


@dataclass(frozen=True)
class MacdaRecord:
    """One selected MACDA time step exactly as it appears in the input file."""

    time_index: int
    mars_date: str
    time_sols: float


@dataclass(frozen=True)
class LabelledRecord:
    """A MACDA record together with the label its output file will carry."""

    time_index: int
    mars_date: MarsDate
    label: WrfPlanetDate
    time_sols: float

    @property
    def hdate(self) -> str:
        return self.label.hdate

    @property
    def filename_stamp(self) -> str:
        return self.label.filename_stamp

    @property
    def relabelled(self) -> bool:
        """True when the output label is not the record's own MACDA date."""

        return self.label != self.mars_date.wrf_planet_date


def parse_mars_date(value: str) -> MarsDate:
    """Parse MACDA ``+00YY-MM-DDThh:mm:ssA`` calendar strings."""

    match = _MARS_DATE_RE.fullmatch(value.strip())
    if match is None:
        raise ValueError(f"Unexpected MACDA Mars_date format: {value!r}")
    if match.group("sign") == "-":
        raise ValueError("Negative Martian years are not supported by MarsWRF HDATE")

    result = MarsDate(year=int(match.group("year")), month=int(match.group("month")),
                      sol=int(match.group("sol")), hour=int(match.group("hour")),
                      minute=int(match.group("minute")),
                      second=int(match.group("second")))
    if not 1 <= result.month <= 12:
        raise ValueError(f"Invalid Martian month in {value!r}")
    max_sol = mars_month_lengths(result.year)[result.month - 1]
    if not 1 <= result.sol <= max_sol:
        raise ValueError(f"Invalid sol {result.sol} for MY {result.year} month "
                         f"{result.month}; valid range is 1..{max_sol}")
    if not 0 <= result.hour <= 23:
        raise ValueError(f"Invalid Martian hour in {value!r}")
    if not 0 <= result.minute <= 59 or not 0 <= result.second <= 59:
        raise ValueError(f"Invalid Martian minute/second in {value!r}")
    return result


def sols_since_macda_epoch(mars_date: MarsDate) -> float:
    """Return sols since MY 24 sol 1 at 00:00 MTC, matching MACDA ``time``."""

    return whole_sols_since_macda_epoch(mars_date) + mars_date.sol_fraction


def whole_sols_since_macda_epoch(mars_date: MarsDate) -> int:
    """Return whole sols since MY 24 sol 1, ignoring the time of sol."""

    return whole_sols_since_epoch(mars_date.year, mars_date.sol_of_year)


def whole_sols_since_epoch(year: int, sol_of_year: int) -> int:
    """Return whole sols from MY 24 sol 1 to ``year`` sol ``sol_of_year``."""

    whole_sols = 0
    if year >= MACDA_EPOCH_YEAR:
        for cycle_year in range(MACDA_EPOCH_YEAR, year):
            whole_sols += mars_year_length(cycle_year)
    else:
        for cycle_year in range(year, MACDA_EPOCH_YEAR):
            whole_sols -= mars_year_length(cycle_year)
    return whole_sols + sol_of_year - 1


def mars_date_at(whole_sols: int, hour: int = 0, minute: int = 0,
                 second: int = 0) -> MarsDate:
    """Inverse of :func:`whole_sols_since_epoch` for MY 24 and later.

    This closes the MACDA calendar loop, so a sol count can be turned back into
    a ``Mars_date`` string without reading a MACDA file.
    """

    if whole_sols < 0:
        raise ValueError(f"Sols before MY {MACDA_EPOCH_YEAR} sol 1 are not supported, "
                         f"got {whole_sols}")
    if not 0 <= hour <= 23:
        raise ValueError(f"Martian hour must be 0..23, got {hour}")
    if not 0 <= minute <= 59 or not 0 <= second <= 59:
        raise ValueError(f"Martian minute/second must be 0..59, got {minute}:{second}")

    year = MACDA_EPOCH_YEAR
    remaining = int(whole_sols)
    while remaining >= mars_year_length(year):
        remaining -= mars_year_length(year)
        year += 1
    month, sol = month_and_sol(year, remaining + 1)
    return MarsDate(year=year, month=month, sol=sol, hour=hour, minute=minute,
                    second=second)


def fixed_year_sols(label: WrfPlanetDate) -> int:
    """Whole sols of a label on the fixed 669-sol calendar metgrid steps through."""

    return label.year * MARSWRF_FIXED_YEAR_SOLS + label.sol_of_year - 1


def label_sols(label: WrfPlanetDate, strategy: str) -> float:
    """Return the continuous sol coordinate of an output label.

    The coordinate is measured on whichever calendar the label is written for:
    ``marswrf_fixed669`` labels live on the fixed 669-sol calendar, every other
    strategy keeps the MACDA sol calendar. Sequence checks compare this value
    against the MACDA ``time`` coordinate, so both sides must use sols.
    """

    check_hdate_strategy(strategy)
    if strategy in ANCHORED_HDATE_STRATEGIES:
        whole_sols = fixed_year_sols(label)
    else:
        whole_sols = whole_sols_since_epoch(label.year, label.sol_of_year)
    return whole_sols + label.sol_fraction


def relabel_fixed_year(mars_date: MarsDate, anchor: MarsDate) -> WrfPlanetDate:
    """Re-label a MACDA date onto the fixed 669-sol MarsWRF calendar.

    ``anchor`` keeps its own MACDA label, and every later sol advances by one
    sol of the fixed 669-sol year. Sol spacing therefore matches what the
    patched WPS ``geth_newdate``/``geth_idts`` compute, so metgrid asks for
    exactly the files a run across a 668-sol MACDA year boundary produces. The
    price is a label offset after the boundary: the MACDA year that WPS thinks
    is one sol longer gains a phantom sol 669, so following records carry a
    label one sol later than their true MACDA date.

    This is a compatibility relabelling only. It does not correct either
    calendar, and it does not move the data in time.
    """

    delta_sols = (whole_sols_since_macda_epoch(mars_date)
                  - whole_sols_since_macda_epoch(anchor))
    anchor_sols = anchor.year * MARSWRF_FIXED_YEAR_SOLS + anchor.sol_of_year - 1
    index = anchor_sols + delta_sols
    if index < 0:
        raise ValueError(f"Fixed {MARSWRF_FIXED_YEAR_SOLS}-sol label for "
                         f"{mars_date.wrf_date} anchored at {anchor.wrf_date} falls "
                         f"before MY 0 sol 1")
    year, sol_index = divmod(index, MARSWRF_FIXED_YEAR_SOLS)
    return WrfPlanetDate(year=year, sol_of_year=sol_index + 1, hour=mars_date.hour,
                         minute=mars_date.minute, second=mars_date.second)


def validate_macda_time(mars_date: MarsDate, time_value: float,
                        tolerance_seconds: float = 1.0) -> None:
    """Ensure independent MACDA ``time`` and ``Mars_date`` coordinates agree."""

    expected = sols_since_macda_epoch(mars_date)
    error_seconds = abs(float(time_value) - expected) * MARS_SECONDS_PER_SOL
    if error_seconds > tolerance_seconds:
        raise ValueError(f"MACDA time mismatch for {mars_date.wrf_date}: "
                         f"time={time_value:.12g} sol, calendar={expected:.12g} sol, "
                         f"error={error_seconds:.3f} s")


def make_label(mars_date: MarsDate, strategy: str = "marswrf_sol",
               anchor: MarsDate | None = None) -> WrfPlanetDate:
    """Return the output label of one record under ``strategy``.

    ``anchor`` is the MACDA date of the run's anchor record, which is the record
    at the configured ``start_index``. ``marswrf_fixed669`` refuses to work
    without it: silently falling back to the record itself would make the same
    MACDA record produce a different HDATE and a different filename depending on
    which time indices one command happened to select.
    """

    check_hdate_strategy(strategy)
    if strategy not in ANCHORED_HDATE_STRATEGIES:
        return mars_date.wrf_planet_date
    if anchor is None:
        raise ValueError(f"hdate_strategy {strategy!r} needs an explicit anchor date; "
                         f"pass the Mars_date of the record at the configured "
                         f"start_index so partial and full runs label alike")
    return relabel_fixed_year(mars_date, anchor)


def make_hdate(mars_date: str, strategy: str = "marswrf_sol",
               anchor_date: str | None = None) -> tuple[str, str]:
    """Return MarsWRF HDATE and the timestamp used in output filenames.

    ``anchor_date`` is the MACDA ``Mars_date`` of the run's anchor record. Only
    ``marswrf_fixed669`` uses it, and that strategy requires it.
    """

    anchor = parse_mars_date(anchor_date) if anchor_date is not None else None
    label = make_label(parse_mars_date(mars_date), strategy, anchor=anchor)
    return label.hdate, label.filename_stamp


def build_label_sequence(records: Sequence[MacdaRecord], strategy: str = "marswrf_sol",
                         anchor_date: str | None = None,
                         tolerance_seconds: float = 1.0) -> list[LabelledRecord]:
    """Label every selected record and validate the sequence as a whole.

    ``records`` are the selected MACDA time steps in write order. The checks are
    the ones a partial rerun needs: labels must be unique so that no output file
    silently overwrites another, strictly increasing so that metgrid reads a
    forward time axis, and spaced exactly like the MACDA ``time`` coordinate so
    that metgrid's own date arithmetic lands on filenames that exist.
    """

    check_hdate_strategy(strategy)
    if not records:
        raise ValueError("No MACDA records selected for labelling")

    anchor = parse_mars_date(anchor_date) if anchor_date is not None else None
    labelled: list[LabelledRecord] = []
    for record in records:
        parsed = parse_mars_date(record.mars_date)
        label = make_label(parsed, strategy, anchor=anchor)
        labelled.append(LabelledRecord(time_index=record.time_index, mars_date=parsed,
                                       label=label, time_sols=float(record.time_sols)))
    check_unique_labels(labelled)
    check_label_spacing(labelled, strategy, tolerance_seconds=tolerance_seconds)
    return labelled


def check_unique_labels(labelled: Sequence[LabelledRecord]) -> None:
    """Reject two records that would be written to the same output filename."""

    seen: dict[str, LabelledRecord] = {}
    for record in labelled:
        stamp = record.filename_stamp
        clash = seen.get(stamp)
        if clash is not None:
            raise ValueError(f"time indices {clash.time_index} and {record.time_index} "
                             f"both label as {stamp}; the second output file would "
                             f"overwrite the first")
        seen[stamp] = record


def check_label_spacing(labelled: Sequence[LabelledRecord], strategy: str,
                        tolerance_seconds: float = 1.0) -> None:
    """Reject labels that are not strictly increasing at the input cadence."""

    tolerance_sols = max(float(tolerance_seconds), 0.0) / MARS_SECONDS_PER_SOL
    previous = labelled[0]
    previous_sols = label_sols(previous.label, strategy)
    for record in labelled[1:]:
        current_sols = label_sols(record.label, strategy)
        if current_sols <= previous_sols:
            raise ValueError(f"labels are not strictly increasing: time index "
                             f"{previous.time_index} labels as "
                             f"{previous.filename_stamp} and {record.time_index} "
                             f"labels as {record.filename_stamp}")
        label_step = current_sols - previous_sols
        time_step = record.time_sols - previous.time_sols
        if abs(label_step - time_step) > tolerance_sols:
            raise ValueError(f"label spacing does not match the MACDA time axis "
                             f"between time indices {previous.time_index} and "
                             f"{record.time_index}: labels {previous.filename_stamp} "
                             f"-> {record.filename_stamp} step {label_step:.12g} sol, "
                             f"time step {time_step:.12g} sol")
        previous = record
        previous_sols = current_sols


def short_mars_years_left(labelled: Sequence[LabelledRecord]) -> list[int]:
    """Return the MACDA years the sequence leaves that are shorter than 669 sols.

    Leaving such a year is exactly the case a fixed 669-sol metgrid cannot
    follow with true MACDA labels: it steps to a sol the year does not have.
    """

    short_years: list[int] = []
    for previous, record in zip(labelled, labelled[1:]):
        year = previous.mars_date.year
        if record.mars_date.year == year:
            continue
        if mars_year_length(year) != MARSWRF_FIXED_YEAR_SOLS and year not in short_years:
            short_years.append(year)
    return short_years
