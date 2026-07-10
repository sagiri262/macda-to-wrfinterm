"""MACDA sol-calendar conversion to the MarsWRF planetary date format."""

from __future__ import annotations

from dataclasses import dataclass
import re


MACDA_EPOCH_YEAR = 24
MARS_SECONDS_PER_SOL = 24 * 60 * 60

# Montabone et al. sol calendar. Month 12 has one extra sol in a 669-sol year.
MONTH_LENGTHS_668 = (56, 55, 56, 55, 56, 56, 55, 56, 55, 56, 56, 56)
_MARS_DATE_RE = re.compile(
    r"^(?P<sign>[+-])(?P<year>\d{4})-(?P<month>\d{2})-(?P<sol>\d{2})"
    r"T(?P<hour>\d{2}):(?P<minute>\d{2}):(?P<second>\d{2})A$"
)


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
    def wrf_date(self) -> str:
        """The 19-character MarsWRF date: YYYY-DDDDD_HH:MM:SS."""

        if not 0 <= self.year <= 9999:
            raise ValueError(f"MarsWRF four-digit year cannot represent MY {self.year}")
        return (
            f"{self.year:04d}-{self.sol_of_year:05d}_"
            f"{self.hour:02d}:{self.minute:02d}:{self.second:02d}"
        )

    @property
    def hdate(self) -> str:
        """The 24-character HDATE stored in WRF intermediate records."""

        return f"{self.wrf_date}.0000"

    @property
    def filename_stamp(self) -> str:
        return self.wrf_date[:13]


def parse_mars_date(value: str) -> MarsDate:
    """Parse MACDA ``+00YY-MM-DDThh:mm:ssA`` calendar strings."""

    match = _MARS_DATE_RE.fullmatch(value.strip())
    if match is None:
        raise ValueError(f"Unexpected MACDA Mars_date format: {value!r}")
    if match.group("sign") == "-":
        raise ValueError("Negative Martian years are not supported by MarsWRF HDATE")

    result = MarsDate(
        year=int(match.group("year")),
        month=int(match.group("month")),
        sol=int(match.group("sol")),
        hour=int(match.group("hour")),
        minute=int(match.group("minute")),
        second=int(match.group("second")),
    )
    if not 1 <= result.month <= 12:
        raise ValueError(f"Invalid Martian month in {value!r}")
    max_sol = mars_month_lengths(result.year)[result.month - 1]
    if not 1 <= result.sol <= max_sol:
        raise ValueError(
            f"Invalid sol {result.sol} for MY {result.year} month {result.month}; "
            f"valid range is 1..{max_sol}"
        )
    if not 0 <= result.hour <= 23:
        raise ValueError(f"Invalid Martian hour in {value!r}")
    if not 0 <= result.minute <= 59 or not 0 <= result.second <= 59:
        raise ValueError(f"Invalid Martian minute/second in {value!r}")
    return result


def sols_since_macda_epoch(mars_date: MarsDate) -> float:
    """Return sols since MY 24 sol 1 at 00:00 MTC, matching MACDA ``time``."""

    whole_sols = 0
    if mars_date.year >= MACDA_EPOCH_YEAR:
        for year in range(MACDA_EPOCH_YEAR, mars_date.year):
            whole_sols += mars_year_length(year)
    else:
        for year in range(mars_date.year, MACDA_EPOCH_YEAR):
            whole_sols -= mars_year_length(year)
    whole_sols += mars_date.sol_of_year - 1
    return whole_sols + mars_date.sol_fraction


def validate_macda_time(
    mars_date: MarsDate,
    time_value: float,
    tolerance_seconds: float = 1.0,
) -> None:
    """Ensure independent MACDA ``time`` and ``Mars_date`` coordinates agree."""

    expected = sols_since_macda_epoch(mars_date)
    error_seconds = abs(float(time_value) - expected) * MARS_SECONDS_PER_SOL
    if error_seconds > tolerance_seconds:
        raise ValueError(
            f"MACDA time mismatch for {mars_date.wrf_date}: time={time_value:.12g} "
            f"sol, calendar={expected:.12g} sol, error={error_seconds:.3f} s"
        )


def make_hdate(mars_date: str, strategy: str = "marswrf_sol") -> tuple[str, str]:
    """Return MarsWRF HDATE and the timestamp used in output filenames."""

    if strategy not in {"marswrf_sol", "mars_date"}:
        raise ValueError(
            f"Unsupported hdate_strategy {strategy!r}; use 'marswrf_sol'"
        )
    parsed = parse_mars_date(mars_date)
    return parsed.hdate, parsed.filename_stamp
