"""MCD season/solar-time sampling on the existing MACDA model calendar.

The model calendar has 24 labelled hours per sol. Its 86400 label seconds
are not SI seconds: physical elapsed time uses 88775.245 SI seconds/sol.
No UTC epoch is inferred. MCD orbit and equation-of-time expressions follow
``MCD_6.2/mcd/MCD.F90`` routines ``sol2ls``, ``ls2sol`` and ``mars_ltime``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterator, Sequence

from .mars_time import MarsDate, mars_month_lengths, mars_year_length


SI_SECONDS_PER_SOL = 88775.245
MODEL_LABEL_SECONDS_PER_SOL = 86400
MCD_ORBIT_SOLS = 668.6
MCD_PERIHELION_SOL = 485.35
MCD_ECCENTRICITY = 0.09340
MCD_PERIHELION_ANGLE = 1.90258341759902
MCD_PERIHELION_LS = 250.99
MCD_OBLIQUITY_DEGREES = 25.1919
MCD_DIURNAL_HOURS = tuple(range(2, 25, 2))


@dataclass(frozen=True)
class SeasonalWeights:
    """
    Monthly-centre interpolation, with zero-based source month indices.

    A nonzero year offset identifies an adjacent-year source that would be
    needed for actual interannual interpolation. If unavailable, a caller
    must explicitly choose climatological cycling or another boundary rule.
    """

    lower_index: int
    upper_index: int
    upper_weight: float
    crosses_year_boundary: bool
    lower_year_offset: int
    upper_year_offset: int

    @property
    def lower_weight(self) -> float:
        return 1.0 - self.upper_weight


@dataclass(frozen=True)
class CyclicWeights:
    """Interpolation into the original order of cyclic diurnal samples."""

    lower_index: int
    upper_index: int
    upper_weight: float

    @property
    def lower_weight(self) -> float:
        return 1.0 - self.upper_weight


@dataclass(frozen=True)
class McdTimePoint:
    """One simultaneous global MTC/LMST0 sample, independent of longitude.

    ``elapsed_sols`` is zero at this model year's sol 1, 00:00. The official
    idealised orbit is 668.6 sols whereas MY36's MACDA label year is 669 sols.
    ``orbit_wrapped`` flags the terminal samples past that orbital period;
    they are not silently rescaled into a longer orbit or mapped to UTC.
    """

    date: MarsDate
    elapsed_sols: float
    elapsed_seconds: float
    calendar_seconds: int
    solar_longitude: float
    mean_solar_time: float
    true_solar_time: float
    seasonal_weights: SeasonalWeights
    orbit_wrapped: bool

    @property
    def wrf_date(self) -> str:
        return self.date.wrf_date

    @property
    def hdate(self) -> str:
        return self.date.hdate

    @property
    def sol_of_year(self) -> int:
        return self.date.sol_of_year


def _finite_number(value: float, name: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _integer(value: int, name: str) -> int:
    result = int(value)
    if result != value or isinstance(value, bool):
        raise ValueError(f"{name} must be an integer")
    return result


def mars_date_from_sol(year: int, sol_of_year: int,
    hour: int = 0, minute: int = 0, second: int = 0) -> MarsDate:
    """Construct the existing MACDA date type from a one-based year sol."""

    # 时间格式
    year = _integer(year, "year")
    sol_of_year = _integer(sol_of_year, "sol_of_year")
    hour = _integer(hour, "hour")
    minute = _integer(minute, "minute")
    second = _integer(second, "second")
    if not 0 <= year <= 9999:
        raise ValueError("year must be representable by four MarsWRF digits")
    if not 1 <= sol_of_year <= mars_year_length(year):
        raise ValueError(f"sol_of_year is outside MY {year}")
    if not 0 <= hour < 24 or not 0 <= minute < 60 or not 0 <= second < 60:
        raise ValueError("hour, minute and second are outside their valid ranges")
    remaining = sol_of_year
    for month, month_length in enumerate(mars_month_lengths(year), start=1):
        if remaining <= month_length:
            return MarsDate(year, month, remaining, hour, minute, second)
        remaining -= month_length
    raise ValueError("calendar month lengths do not span the requested sol")


def sol_to_ls(sol: float) -> float:
    """Return MCD idealised Ls [degrees] from sols after model spring equinox.

    This follows the official Kepler solution in double precision. Its
    rounded orbital constants give a small equinox offset (about 0.000063 deg),
    which is retained instead of forcing Ls=0 for the first timestamp.
    """

    sol = _finite_number(sol, "sol")
    phase = (sol - MCD_PERIHELION_SOL) / MCD_ORBIT_SOLS
    mean_anomaly = 2.0 * math.pi * (phase - math.floor(phase + 0.5))
    reference = abs(mean_anomaly)
    eccentric_anomaly = reference + MCD_ECCENTRICITY * math.sin(reference)
    for iteration in range(10):
        residual = eccentric_anomaly - MCD_ECCENTRICITY * math.sin(eccentric_anomaly) - reference
        derivative = 1.0 - MCD_ECCENTRICITY * math.cos(eccentric_anomaly)
        correction = -residual / derivative
        eccentric_anomaly += correction
        if abs(correction) <= 1.0e-7:
            break
    if mean_anomaly < 0:
        eccentric_anomaly = -eccentric_anomaly
    orbit_factor = math.sqrt((1.0 + MCD_ECCENTRICITY) / (1.0 - MCD_ECCENTRICITY))
    true_anomaly = 2.0 * math.atan(orbit_factor * math.tan(eccentric_anomaly / 2.0))
    return math.degrees(true_anomaly - MCD_PERIHELION_ANGLE) % 360.0


def ls_to_sol(solar_longitude: float) -> float:
    """Return MCD model sol for Ls [degrees], following official ``ls2sol``."""

    solar_longitude = _finite_number(solar_longitude, "solar_longitude")
    if not 0.0 <= solar_longitude <= 360.0:
        raise ValueError("solar_longitude must be between 0 and 360 degrees")
    if abs(solar_longitude) < 1.0e-5:
        return 0.0
    true_anomaly = math.radians(solar_longitude) + MCD_PERIHELION_ANGLE
    orbit_factor = math.sqrt((1.0 + MCD_ECCENTRICITY) / (1.0 - MCD_ECCENTRICITY))
    eccentric_anomaly = 2.0 * math.atan(math.tan(true_anomaly / 2.0) / orbit_factor)
    mean_anomaly = eccentric_anomaly - MCD_ECCENTRICITY * math.sin(eccentric_anomaly)
    model_sol = mean_anomaly * MCD_ORBIT_SOLS / (2.0 * math.pi) + MCD_PERIHELION_SOL
    return model_sol % MCD_ORBIT_SOLS


def equation_of_time_hours(solar_longitude: float) -> float:
    """Return MCD E_T in Martian hours, with LTST0 = LMST0 - E_T."""

    solar_longitude = _finite_number(solar_longitude, "solar_longitude")
    orbital_term = 2.0 * MCD_ECCENTRICITY * math.sin(math.radians(solar_longitude - MCD_PERIHELION_LS))
    axial_factor = math.tan(math.radians(MCD_OBLIQUITY_DEGREES) / 2.0) ** 2
    axial_term = axial_factor * math.sin(2.0 * math.radians(solar_longitude))
    return (orbital_term - axial_term) * 24.0 / (2.0 * math.pi)


def mean_to_true_solar_time(mean_solar_hour: float,
    solar_longitude: float, longitude: float = 0.0) -> float:
    """Convert prime-meridian LMST to LTST at east-positive longitude.

    Use longitude=0 when querying an entire MCD map for one physical time.
    Supplying each grid longitude would produce local times, not a common
    prime-meridian sampling index for the global MCD array.
    """

    mean_solar_hour = _finite_number(mean_solar_hour, "mean_solar_hour")
    longitude = _finite_number(longitude, "longitude")
    return (mean_solar_hour - equation_of_time_hours(solar_longitude) + longitude / 15.0) % 24.0


def monthly_weights(solar_longitude: float) -> SeasonalWeights:
    """Interpolate between 15,45,...,345-degree Ls month centres."""

    solar_longitude = _finite_number(solar_longitude, "solar_longitude") % 360.0
    centred_phase = (solar_longitude - 15.0) / 30.0
    lower_unwrapped = math.floor(centred_phase)
    upper_unwrapped = lower_unwrapped + 1
    upper_weight = centred_phase - lower_unwrapped
    lower_index = lower_unwrapped % 12
    upper_index = upper_unwrapped % 12
    lower_offset = lower_unwrapped // 12
    upper_offset = upper_unwrapped // 12
    crosses_boundary = lower_offset != upper_offset
    return SeasonalWeights(lower_index, upper_index, upper_weight,
        crosses_boundary, lower_offset, upper_offset)


def diurnal_weights(solar_hour: float,
    sample_hours: Sequence[float] = MCD_DIURNAL_HOURS) -> CyclicWeights:
    """Compute cyclic linear weights while preserving source array order.

    MCD normally supplies 2,4,...,24 h. The final 24 h sample is the 0 h
    endpoint; 23 h therefore uses 22/24 h, and 1 h uses 24/2 h.
    """

    solar_hour = _finite_number(solar_hour, "solar_hour") % 24.0
    if len(sample_hours) < 2:
        raise ValueError("at least two diurnal samples are required")
    samples = [(_finite_number(value, "sample_hours") % 24.0, index)
        for index, value in enumerate(sample_hours)]
    samples.sort()
    if any(samples[index][0] == samples[index - 1][0] for index in range(1, len(samples))):
        raise ValueError("diurnal samples must be unique modulo 24 hours")
    for index, (lower_hour, lower_index) in enumerate(samples):
        upper_hour, upper_index = samples[(index + 1) % len(samples)]
        if index == len(samples) - 1:
            upper_hour += 24.0
        query_hour = solar_hour
        if query_hour < lower_hour:
            query_hour += 24.0
        if lower_hour <= query_hour < upper_hour:
            upper_weight = (query_hour - lower_hour) / (upper_hour - lower_hour)
            return CyclicWeights(lower_index, upper_index, upper_weight)
    raise ValueError("failed to locate solar hour in the cyclic sample grid")


def iter_mcd_times(year: int = 36, step_hours: int = 2,
    start_sol: int = 1, end_sol: int | None = None) -> Iterator[McdTimePoint]:
    """Yield an inclusive sol range, sampled at labelled prime-meridian LMST.

    MY36 defaults to 669 sols and 8028 samples from 00:00 through 22:00.
    Every sample is a simultaneous global state. The calendar's labels are
    not UTC; elapsed_seconds is a relative physical SI-second coordinate.
    """

    year = _integer(year, "year")
    step_hours = _integer(step_hours, "step_hours")
    start_sol = _integer(start_sol, "start_sol")
    end_sol = mars_year_length(year) if end_sol is None else _integer(end_sol, "end_sol")
    if not 0 <= year <= 9999:
        raise ValueError("year must be representable by four MarsWRF digits")
    if step_hours <= 0 or 24 % step_hours != 0:
        raise ValueError("step_hours must be a positive integer divisor of 24")
    if not 1 <= start_sol <= end_sol <= mars_year_length(year):
        raise ValueError("requested sol range is outside the selected calendar year")
    for sol_of_year in range(start_sol, end_sol + 1):
        for hour in range(0, 24, step_hours):
            elapsed_sols = sol_of_year - 1 + hour / 24.0
            solar_longitude = sol_to_ls(elapsed_sols)
            mean_solar_hour = float(hour)
            true_solar_hour = mean_to_true_solar_time(mean_solar_hour, solar_longitude)
            date = mars_date_from_sol(year, sol_of_year, hour)
            calendar_seconds = (sol_of_year - 1) * MODEL_LABEL_SECONDS_PER_SOL + hour * 3600
            yield McdTimePoint(date, elapsed_sols, elapsed_sols * SI_SECONDS_PER_SOL,
                calendar_seconds, solar_longitude, mean_solar_hour, true_solar_hour,
                monthly_weights(solar_longitude), elapsed_sols >= MCD_ORBIT_SOLS)
