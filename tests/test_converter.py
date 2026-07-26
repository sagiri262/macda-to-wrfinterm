from __future__ import annotations

from types import SimpleNamespace
import unittest

from macda2wrf.converter import MacdaConverter
from tests.macda_boundary_fixture import (FIRST_MY35_INDEX, LAST_INDEX, LAST_MY34_INDEX,
                                          RECORD_COUNT, FakeMacdaReader,
                                          boundary_records, time_axis)


def make_converter(strategy: str = "marswrf_fixed669", start_index: int = 0,
                   end_index: int | None = None, max_times: int | None = None,
                   xfcst_from_start: bool = True) -> MacdaConverter:
    """Build a converter with only the config a labelling test touches.

    ``object.__new__`` skips ``__init__`` so no MACDA file, variable table, or
    target grid is needed.
    """

    converter = object.__new__(MacdaConverter)
    converter.cfg = SimpleNamespace(start_index=start_index, end_index=end_index,
                                    max_times=max_times, hdate_strategy=strategy,
                                    time_tolerance_seconds=1.0,
                                    xfcst_from_start=xfcst_from_start)
    return converter


class ConverterIndexTest(unittest.TestCase):
    def setUp(self) -> None:
        self.converter = make_converter(max_times=1)

    def test_config_limit(self) -> None:
        self.assertEqual(self.converter._time_indices(4), [0])

    def test_cli_override_can_expand_config_limit(self) -> None:
        self.assertEqual(self.converter._time_indices(4, apply_config_limit=False),
                         [0, 1, 2, 3])

    def test_explicit_indices_are_sorted_and_deduplicated(self) -> None:
        converter = make_converter()
        self.assertEqual(converter._select_indices(RECORD_COUNT, 0, [47, 0, 47], None),
                         [0, 47])

    def test_indices_before_the_anchor_are_rejected(self) -> None:
        converter = make_converter(start_index=10)
        with self.assertRaisesRegex(ValueError, "come before the anchor record"):
            converter._select_indices(RECORD_COUNT, 10, [5, 20], None)

    def test_indices_outside_the_file_are_rejected(self) -> None:
        converter = make_converter()
        with self.assertRaisesRegex(ValueError, "outside available range"):
            converter._select_indices(RECORD_COUNT, 0, [0, RECORD_COUNT], None)

    def test_anchor_must_exist_in_the_file(self) -> None:
        converter = make_converter(start_index=RECORD_COUNT)
        with self.assertRaisesRegex(ValueError, "label and XFCST anchor"):
            converter._anchor_index(RECORD_COUNT)


class ConverterStableAnchorTest(unittest.TestCase):
    """A partial rerun must write the filenames and XFCST of the full run."""

    def setUp(self) -> None:
        self.records = boundary_records()
        self.reader = FakeMacdaReader(self.records)
        self.time_values = time_axis(self.records)

    def labels(self, converter: MacdaConverter, selection) -> list:
        return converter._build_labels(self.reader, list(selection), self.time_values, 0)

    def test_single_index_rerun_matches_the_full_run(self) -> None:
        converter = make_converter()
        full = {record.time_index: record.hdate
                for record in self.labels(converter, range(RECORD_COUNT))}
        for index in (0, LAST_MY34_INDEX, FIRST_MY35_INDEX, 58, 59, LAST_INDEX):
            single = self.labels(converter, [index])
            self.assertEqual(single[0].hdate, full[index])

    def test_index_47_alone_still_gets_the_phantom_sol(self) -> None:
        converter = make_converter()
        single = self.labels(converter, [FIRST_MY35_INDEX])
        self.assertEqual(single[0].filename_stamp, "0034-00669_00")
        self.assertEqual(single[0].mars_date.filename_stamp, "0035-00001_00")
        self.assertTrue(single[0].relabelled)

    def test_default_strategy_keeps_the_true_label(self) -> None:
        converter = make_converter(strategy="marswrf_sol")
        single = self.labels(converter, [FIRST_MY35_INDEX])
        self.assertEqual(single[0].filename_stamp, "0035-00001_00")
        self.assertFalse(single[0].relabelled)

    def test_anchor_is_the_config_start_index_not_the_selection(self) -> None:
        converter = make_converter(start_index=0)
        # A rerun that only asks for the two records around the boundary must
        # still place them on the fixed calendar of the full run.
        pair = self.labels(converter, [LAST_MY34_INDEX, FIRST_MY35_INDEX])
        self.assertEqual([record.filename_stamp for record in pair],
                         ["0034-00668_22", "0034-00669_00"])

    def test_xfcst_is_relative_to_the_config_anchor(self) -> None:
        converter = make_converter()
        for index in (0, FIRST_MY35_INDEX, LAST_INDEX):
            hours = converter._xfcst_hours(self.time_values, index, 0)
            self.assertAlmostEqual(hours, 2.0 * index)

    def test_xfcst_does_not_depend_on_the_selection(self) -> None:
        converter = make_converter()
        self.assertAlmostEqual(converter._xfcst_hours(self.time_values, 47, 0), 94.0)

    def test_xfcst_before_the_anchor_is_rejected(self) -> None:
        converter = make_converter(start_index=FIRST_MY35_INDEX)
        with self.assertRaisesRegex(ValueError, "before the anchor record"):
            converter._xfcst_hours(self.time_values, 0, FIRST_MY35_INDEX)

    def test_xfcst_can_be_switched_off(self) -> None:
        converter = make_converter(xfcst_from_start=False)
        self.assertEqual(converter._xfcst_hours(self.time_values, 47, 0), 0.0)


if __name__ == "__main__":
    unittest.main()
