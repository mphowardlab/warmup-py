import warmup
import numpy as np
import unittest

class TestRadialDistribution(unittest.TestCase):
    """Tests for the RadialDistribution function calculator."""

    def test_init(self):
        """Test RadialDistribution initialization."""
        rdf = warmup.RadialDistribution(rmax=1.0, dr=0.1)
        self.assertEqual(rdf.rmin, 0.0)
        self.assertEqual(rdf.rmax, 1.0)
        self.assertEqual(rdf.dr, 0.1)

        np.testing.assert_allclose(rdf.edges, np.arange(0.0,1.05,0.1))
        np.testing.assert_allclose(rdf.centers, np.arange(0.05,1.05,0.1))
        np.testing.assert_allclose(rdf.rdf, np.zeros(10))

    def test_range(self):
        """Test binning range for RadialDistribution."""
        rdf = warmup.RadialDistribution(rmin=0.5, rmax=1.5, dr=0.2)
        self.assertEqual(rdf.rmin, 0.5)
        self.assertEqual(rdf.rmax, 1.5)
        self.assertEqual(rdf.dr, 0.2)

        np.testing.assert_allclose(rdf.edges, np.arange(0.5,1.6,0.2))
        np.testing.assert_allclose(rdf.centers, np.arange(0.6,1.6,0.2))
        np.testing.assert_allclose(rdf.rdf, np.zeros(5))

    def test_nofit_dr(self):
        """Test for bin width that does not fit."""
        rdf = warmup.RadialDistribution(rmax=1.0,dr=0.12)
        self.assertAlmostEqual(rdf.dr, 0.125)

        rdf = warmup.RadialDistribution(rmax=1.0,dr=0.101)
        self.assertAlmostEqual(rdf.dr, 0.1)

    def test_accumulate(self):
        """Test accumulation for RadialDistribution on simple cubic lattice."""
        # simple cubic lattice
        L = 5.0
        a = 1.0
        n = np.round(L/a).astype(int)
        points = [[a*i,a*j,a*k] for i in range(n) for j in range(n) for k in range(n)]

        # only count the first nearest neighbor
        rdf = warmup.RadialDistribution(rmax=1.01, dr=0.101)
        gsc = np.zeros(len(rdf.centers))
        gsc[-1] = 5.13014 # from ovito / analytical formula

        # run once, should give target
        rdf.accumulate(L,points)
        np.testing.assert_allclose(rdf.rdf, gsc, rtol=1.e-5)

        # run a second time, should give same result
        rdf.accumulate(L,points)
        np.testing.assert_allclose(rdf.rdf, gsc, rtol=1.e-5)

        # reset to zero
        rdf.reset()
        np.testing.assert_allclose(rdf.rdf, np.zeros(len(rdf.centers)))

        # rerun a third time
        rdf.accumulate(L,points)
        np.testing.assert_allclose(rdf.rdf, gsc, rtol=1.e-5)

    def test_minimage_error(self):
        """Test RadialDistribution function checks minimum image."""
        rdf = warmup.RadialDistribution(rmax=1.0, dr=0.1)
        with self.assertRaises(ValueError):
            rdf.accumulate(1.9, np.zeros((2,3)))

if __name__ == '__main__':
    unittest.main()
