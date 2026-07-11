from __future__ import annotations

import unittest

from macda2wrf.mars_time import (
    make_hdate,
    mars_month_lengths,
    mars_year_length,
    parse_mars_date,
    sols_since_macda_epoch,
    validate_macda_time,
)


class MarsTimeTest(unittest.TestCase):
    def test_five_year_cycle(self) -> None:
        self.assertEqual(
            [mars_year_length(year) for year in range(21, 26)],
            [669, 668, 669, 668, 669],
        )
        self.assertEqual(sum(mars_month_lengths(24)), 668)
        self.assertEqual(sum(mars_month_lengths(28)), 669)

    def test_requested_file_first_time(self) -> None:
        mars_date = parse_mars_date("+0028-10-07T02:00:00A")
        self.assertEqual(mars_date.sol_of_year, 507)
        self.assertEqual(mars_date.wrf_date, "0028-00507_02:00:00")
        self.assertAlmostEqual(sols_since_macda_epoch(mars_date), 3180.083333333333)
        validate_macda_time(mars_date, 3180.08333333333)

    def test_requested_file_last_time(self) -> None:
        hdate, stamp = make_hdate("+0028-10-37T00:00:00A")
        self.assertEqual(hdate, "0028-00537_00:00:00.0000")
        self.assertEqual(stamp, "0028-00537_00")

    def test_rejects_bad_sol_and_mismatch(self) -> None:
        with self.assertRaisesRegex(ValueError, "Invalid sol"):
            parse_mars_date("+0028-02-56T00:00:00A")
        with self.assertRaisesRegex(ValueError, "time mismatch"):
            validate_macda_time(
                parse_mars_date("+0028-10-07T02:00:00A"),
                3181.0,
            )


if __name__ == "__main__":
    unittest.main()
