from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from macda2wrf.grid import RegularLatLonGrid
from macda2wrf.wrf_intermediate import (
    WrfIntermediateWriter,
    read_intermediate_file,
    validate_intermediate_file,
)


class IntermediateRoundTripTest(unittest.TestCase):
    def test_round_trip(self) -> None:
        grid = RegularLatLonGrid(-2.5, 2.5, 2, 0.0, 5.0, 2, 3389.92)
        slab = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "MACDA:0028-00507_02"
            with WrfIntermediateWriter(
                path,
                grid,
                "0028-00507_02:00:00.0000",
                "MACDA",
            ) as writer:
                writer.write_field("UU", slab, "m s-1", "Zonal wind", 610.0, 1)

            records = read_intermediate_file(path)
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0].field, "UU")
            self.assertEqual(records[0].hdate, "0028-00507_02:00:00.0000")
            self.assertEqual(records[0].xfcst, 0.0)
            np.testing.assert_array_equal(records[0].slab, slab)
            validate_intermediate_file(
                path,
                expected_hdate="0028-00507_02:00:00.0000",
                expected_shape=(2, 2),
                required_fields={"UU"},
            )


if __name__ == "__main__":
    unittest.main()
