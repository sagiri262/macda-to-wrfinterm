from __future__ import annotations

from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from utils.extractNetcdfYear import extractYear


def writeSourceFile(filePath: Path, hours: list[float]) -> None:
    with Dataset(filePath, "w", format="NETCDF4") as dataset:
        dataset.createDimension("time", None)
        dataset.createDimension("latitude", 2)
        timeVariable = dataset.createVariable("time", "f8", ("time",))
        timeVariable.units = "hours since 2000-11-30 00:00:00"
        timeVariable.calendar = "standard"
        temperatureVariable = dataset.createVariable(
            "temperature",
            "f4",
            ("time", "latitude"),
        )
        latitudeVariable = dataset.createVariable("latitude", "f4", ("latitude",))
        timeVariable[:] = hours
        temperatureVariable[:] = np.asarray(
            [[hour, hour + 0.5] for hour in hours],
            dtype=np.float32,
        )
        latitudeVariable[:] = [10.0, 20.0]


def testExtractYearSortsFilesAndDropsBoundaryDuplicates(tmp_path: Path) -> None:
    firstSourcePath = tmp_path / "z_cross_year.nc"
    secondSourcePath = tmp_path / "a_cross_year.nc"
    outputPath = tmp_path / "reanalysis_2001.nc"
    hourPerDay = 24.0
    firstYearHour = 32 * hourPerDay
    writeSourceFile(
        firstSourcePath,
        [firstYearHour - hourPerDay, firstYearHour, firstYearHour + hourPerDay],
    )
    writeSourceFile(
        secondSourcePath,
        [firstYearHour + hourPerDay, firstYearHour + 2 * hourPerDay, 10000.0],
    )

    recordCount, duplicateCount = extractYear(
        [secondSourcePath, firstSourcePath],
        outputPath,
        2001,
        None,
        1,
    )

    assert recordCount == 3
    assert duplicateCount == 1
    with Dataset(outputPath, "r") as dataset:
        assert dataset.variables["time"][:].tolist() == [768.0, 792.0, 816.0]
        assert dataset.variables["temperature"][:, 0].tolist() == [768.0, 792.0, 816.0]
        assert dataset.variables["latitude"][:].tolist() == [10.0, 20.0]
