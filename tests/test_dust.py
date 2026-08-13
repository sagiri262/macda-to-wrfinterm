from __future__ import annotations

import unittest

import numpy as np

from macda2wrf.dust import (columnOpticalDepth, combineDustTracers, opticalDepthScaling,
                            sigmaLayerThickness, tracerFromDustMass, validateDustParameters,
                            validateDustPartition)


class DustMassTest(unittest.TestCase):
    def testPartitionPreservesTotalMass(self) -> None:
        source = np.array([1.0e-6, 3.0e-6], dtype=np.float64)
        tracer1 = tracerFromDustMass(source, 0.2, 0.3)
        tracer2 = tracerFromDustMass(source, 0.8, 0.7)
        validateDustPartition(source, tracer1, tracer2, 0.3, 0.7)
        np.testing.assert_allclose(combineDustTracers(tracer1, tracer2, 0.3, 0.7), source)

    def testInvalidFractionsAreRejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "must sum to 1.0"):
            validateDustParameters(0.2, 0.7, 0.3, 0.7, (1.0,))


class DustOpticsTest(unittest.TestCase):
    def testSigmaThicknessCoversWholeColumn(self) -> None:
        thickness = sigmaLayerThickness(np.array([0.9, 0.5, 0.1]))
        self.assertAlmostEqual(float(np.sum(thickness)), 1.0)
        self.assertTrue(np.all(thickness > 0.0))

    def testUniformDustColumnHasAnalyticOpticalDepth(self) -> None:
        dust = np.full((3, 2, 2), 2.0e-6)
        pressure = np.full((2, 2), 600.0)
        opticalDepth = columnOpticalDepth(dust, pressure, np.array([0.9, 0.5, 0.1]), 500.0, 3.72)
        expected = np.full((2, 2), 500.0 * 600.0 * 2.0e-6 / 3.72)
        np.testing.assert_allclose(opticalDepth, expected)

    def testCurrentMarswrfOpticsDoNotCloseMacdaOpacity(self) -> None:
        scale = opticalDepthScaling(2500.0, 1.5e-6, 2.62718, 3.72, 2650.0, 3.727,
            1.0e-6, 2.4e-6, 3.04, 3.04, 0.3, 0.7)
        self.assertGreater(scale, 1.0)
        self.assertFalse(np.isclose(scale, 1.0, rtol=0.05, atol=0.0))


if __name__ == "__main__":
    unittest.main()
