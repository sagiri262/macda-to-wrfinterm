from __future__ import annotations

from types import SimpleNamespace
import unittest

from macda2wrf.converter import MacdaConverter


class ConverterIndexTest(unittest.TestCase):
    def setUp(self) -> None:
        self.converter = object.__new__(MacdaConverter)
        self.converter.cfg = SimpleNamespace(
            start_index=0,
            end_index=None,
            max_times=1,
        )

    def test_config_limit(self) -> None:
        self.assertEqual(self.converter._time_indices(4), [0])

    def test_cli_override_can_expand_config_limit(self) -> None:
        self.assertEqual(
            self.converter._time_indices(4, apply_config_limit=False),
            [0, 1, 2, 3],
        )


if __name__ == "__main__":
    unittest.main()
