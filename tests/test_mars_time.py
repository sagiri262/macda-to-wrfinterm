from __future__ import annotations

import unittest

from macda2wrf.mars_time import (MacdaRecord, build_label_sequence, label_sols,
                                 make_hdate, mars_date_at, mars_month_lengths,
                                 mars_year_length, parse_mars_date, relabel_fixed_year,
                                 short_mars_years_left, sols_since_macda_epoch,
                                 validate_macda_time, whole_sols_since_epoch)
from tests.macda_boundary_fixture import (FIRST_MY35_INDEX, LAST_INDEX,
                                          LAST_MY34_INDEX, RECORD_COUNT,
                                          boundary_records)


class MarsTimeTest(unittest.TestCase):
    def test_five_year_cycle(self) -> None:
        self.assertEqual([mars_year_length(year) for year in range(21, 26)],
                         [669, 668, 669, 668, 669])
        self.assertEqual(sum(mars_month_lengths(24)), 668)
        self.assertEqual(sum(mars_month_lengths(28)), 669)
        self.assertEqual(mars_year_length(34), 668)
        self.assertEqual(mars_year_length(35), 669)

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
            validate_macda_time(parse_mars_date("+0028-10-07T02:00:00A"), 3181.0)


class MacdaCalendarInverseTest(unittest.TestCase):
    """The synthetic fixture must reproduce the real cross-year file exactly."""

    def test_fixture_endpoints_match_the_real_file(self) -> None:
        records = boundary_records()
        self.assertEqual(len(records), RECORD_COUNT)
        self.assertEqual(records[0].mars_date, "+0034-12-53T02:00:00A")
        self.assertEqual(parse_mars_date(records[0].mars_date).wrf_date,
                         "0034-00665_02:00:00")
        self.assertEqual(parse_mars_date(records[LAST_INDEX].mars_date).wrf_date,
                         "0035-00027_00:00:00")

    def test_fixture_reproduces_the_macda_year_boundary(self) -> None:
        records = boundary_records()
        self.assertEqual(parse_mars_date(records[LAST_MY34_INDEX].mars_date).wrf_date,
                         "0034-00668_22:00:00")
        self.assertEqual(parse_mars_date(records[FIRST_MY35_INDEX].mars_date).wrf_date,
                         "0035-00001_00:00:00")

    def test_mars_date_at_inverts_the_sol_count(self) -> None:
        for record in boundary_records():
            parsed = parse_mars_date(record.mars_date)
            whole_sols = whole_sols_since_epoch(parsed.year, parsed.sol_of_year)
            self.assertEqual(mars_date_at(whole_sols, hour=parsed.hour), parsed)
            self.assertEqual(parsed.macda_string, record.mars_date)

    def test_fixture_time_axis_matches_its_own_calendar(self) -> None:
        for record in boundary_records():
            validate_macda_time(parse_mars_date(record.mars_date), record.time_sols)


class FixedYearLabelTest(unittest.TestCase):
    """MY 34 (668 sols) to MY 35, anchored at MY 34 sol 665 02:00 MTC."""

    ANCHOR = "+0034-12-53T02:00:00A"

    def label(self, mars_date: str) -> str:
        hdate, stamp = make_hdate(mars_date, "marswrf_fixed669", anchor_date=self.ANCHOR)
        self.assertEqual(stamp, hdate[:13])
        return hdate[:19]

    def test_anchor_keeps_its_macda_label(self) -> None:
        self.assertEqual(self.label(self.ANCHOR), "0034-00665_02:00:00")

    def test_labels_stay_true_before_the_year_boundary(self) -> None:
        self.assertEqual(self.label("+0034-12-56T22:00:00A"), "0034-00668_22:00:00")

    def test_boundary_gets_the_phantom_sol_669(self) -> None:
        # True MY 35 sol 1 is labelled MY 34 sol 669, the sol WPS expects but
        # the 668-sol MACDA calendar does not have.
        self.assertEqual(self.label("+0035-01-01T00:00:00A"), "0034-00669_00:00:00")

    def test_labels_after_the_boundary_lag_by_one_sol(self) -> None:
        self.assertEqual(self.label("+0035-01-02T00:00:00A"), "0035-00001_00:00:00")
        self.assertEqual(self.label("+0035-01-27T00:00:00A"), "0035-00026_00:00:00")

    def test_labels_advance_one_fixed_year_sol_per_real_sol(self) -> None:
        anchor = parse_mars_date(self.ANCHOR)
        sequence = ["+0034-12-53T02:00:00A", "+0034-12-56T22:00:00A",
                    "+0035-01-01T00:00:00A", "+0035-01-02T00:00:00A",
                    "+0035-01-27T00:00:00A"]
        for value in sequence:
            parsed = parse_mars_date(value)
            label = relabel_fixed_year(parsed, anchor)
            # This is what the patched WPS geth_idts/geth_newdate compute.
            wps_sols = label.year * 669 + label.sol_of_year
            real_sols = sols_since_macda_epoch(parsed)
            self.assertAlmostEqual(wps_sols + parsed.sol_fraction - real_sols,
                                   669 * 34 + 665 - sols_since_macda_epoch(anchor)
                                   + anchor.sol_fraction)

    def test_single_year_run_is_unchanged(self) -> None:
        hdate, _ = make_hdate("+0028-10-37T00:00:00A", "marswrf_fixed669",
                              anchor_date="+0028-10-07T02:00:00A")
        self.assertEqual(hdate, "0028-00537_00:00:00.0000")

    def test_rejects_unknown_strategy(self) -> None:
        with self.assertRaisesRegex(ValueError, "Unsupported hdate_strategy"):
            make_hdate("+0028-10-07T02:00:00A", "gregorian")

    def test_requires_an_explicit_anchor(self) -> None:
        with self.assertRaisesRegex(ValueError, "needs an explicit anchor date"):
            make_hdate("+0035-01-01T00:00:00A", "marswrf_fixed669")

    def test_default_strategy_keeps_the_true_my35_label(self) -> None:
        hdate, stamp = make_hdate("+0035-01-01T00:00:00A")
        self.assertEqual(hdate, "0035-00001_00:00:00.0000")
        self.assertEqual(stamp, "0035-00001_00")


class StableAnchorLabelTest(unittest.TestCase):
    """One command's selection must never move the labels of the other records."""

    STRATEGY = "marswrf_fixed669"

    def setUp(self) -> None:
        self.records = boundary_records()
        self.anchor_date = self.records[0].mars_date

    def labels(self, selection, strategy: str | None = None) -> list:
        picked = [self.records[index] for index in selection]
        return build_label_sequence(picked, strategy=strategy or self.STRATEGY,
                                    anchor_date=self.anchor_date)

    def test_single_record_rerun_matches_the_full_run(self) -> None:
        full = {record.time_index: record.hdate
                for record in self.labels(range(RECORD_COUNT))}
        for index in (0, LAST_MY34_INDEX, FIRST_MY35_INDEX, 58, 59, LAST_INDEX):
            single = self.labels([index])
            self.assertEqual(single[0].hdate, full[index])
            self.assertEqual(single[0].filename_stamp, full[index][:13])

    def test_full_run_labels_the_year_boundary_on_the_fixed_calendar(self) -> None:
        labels = {record.time_index: record.hdate[:19]
                  for record in self.labels(range(RECORD_COUNT))}
        self.assertEqual(labels[LAST_MY34_INDEX], "0034-00668_22:00:00")
        self.assertEqual(labels[FIRST_MY35_INDEX], "0034-00669_00:00:00")
        self.assertEqual(labels[58], "0034-00669_22:00:00")
        self.assertEqual(labels[59], "0035-00001_00:00:00")
        self.assertEqual(labels[LAST_INDEX], "0035-00026_00:00:00")

    def test_labels_are_continuous_across_the_boundary(self) -> None:
        labelled = self.labels(range(RECORD_COUNT))
        coordinates = [label_sols(record.label, self.STRATEGY) for record in labelled]
        for previous, current in zip(coordinates, coordinates[1:]):
            self.assertAlmostEqual(current - previous, 2.0 / 24.0)

    def test_default_strategy_keeps_every_true_macda_label(self) -> None:
        labelled = self.labels(range(RECORD_COUNT), strategy="marswrf_sol")
        labels = {record.time_index: record.hdate[:19] for record in labelled}
        self.assertEqual(labels[FIRST_MY35_INDEX], "0035-00001_00:00:00")
        self.assertEqual(labels[LAST_INDEX], "0035-00027_00:00:00")
        self.assertFalse(any(record.relabelled for record in labelled))

    def test_short_mars_year_is_reported_for_true_macda_labels(self) -> None:
        labelled = self.labels(range(RECORD_COUNT), strategy="marswrf_sol")
        self.assertEqual(short_mars_years_left(labelled), [34])

    def test_sequence_requires_an_anchor_for_the_fixed_calendar(self) -> None:
        with self.assertRaisesRegex(ValueError, "needs an explicit anchor date"):
            build_label_sequence(self.records[:4], strategy=self.STRATEGY)


class LabelSequenceValidationTest(unittest.TestCase):
    def test_duplicate_labels_are_rejected(self) -> None:
        records = boundary_records(2)
        clash = MacdaRecord(time_index=1, mars_date=records[0].mars_date,
                            time_sols=records[0].time_sols)
        with self.assertRaisesRegex(ValueError, "overwrite the first"):
            build_label_sequence([records[0], clash])

    def test_labels_must_be_strictly_increasing(self) -> None:
        records = boundary_records(3)
        with self.assertRaisesRegex(ValueError, "not strictly increasing"):
            build_label_sequence([records[0], records[2], records[1]])

    def test_label_spacing_must_match_the_time_axis(self) -> None:
        records = boundary_records(2)
        drifted = MacdaRecord(time_index=1, mars_date=records[1].mars_date,
                              time_sols=records[1].time_sols + 0.5)
        with self.assertRaisesRegex(ValueError, "label spacing does not match"):
            build_label_sequence([records[0], drifted])

    def test_empty_selection_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "No MACDA records selected"):
            build_label_sequence([])

    def test_unknown_strategy_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "Unsupported hdate_strategy"):
            build_label_sequence(boundary_records(2), strategy="gregorian")


if __name__ == "__main__":
    unittest.main()
