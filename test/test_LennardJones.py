import warmup
import numpy as np
import unittest

class TestLennardJones(unittest.TestCase):
    """Tests for the Lennard-Jones pair potential."""

    def setUp(self):
        self.lj = warmup.LennardJones(epsilon=1.5, sigma=0.5, rcut=2.5)

    def test_init(self):
        """Test Lennard-Jones constructor."""
        self.assertEqual(self.lj.epsilon, 1.5)
        self.assertEqual(self.lj.sigma, 0.5)
        self.assertEqual(self.lj.rcut, 2.5)

    def test_energy(self):
        """Test Lennard-Jones energy."""
        self.assertAlmostEqual(self.lj.energy(0.5), 0.0)
        self.assertAlmostEqual(self.lj.energy(0.5*2**(1./6.)), -1.5)
        self.assertAlmostEqual(self.lj.energy(3.0), 0.0)

        r = [0.5,1.0,3.0]
        u = self.lj.energy(r)
        np.testing.assert_allclose(u,[0,-63*1.5/1024,0])

    def test_force(self):
        """Test Lennard-Jones force."""
        self.assertAlmostEqual(self.lj.force(0.5), 24*1.5/0.5)
        self.assertAlmostEqual(self.lj.force(0.5*2**(1./6.)), 0.0)
        self.assertAlmostEqual(self.lj.force(3.0), 0.0)

        r = [0.5,1.0,3.0]
        f = self.lj.force(r)
        np.testing.assert_allclose(f,[24*1.5/0.5,-(93./512.)*1.5/0.5,0])

    def test_inf(self):
        """Test infinities are raised at r=0."""
        self.assertEqual(self.lj.energy(0.0), np.inf)
        self.assertEqual(self.lj.force(0.0), np.inf)

if __name__ == '__main__':
    unittest.main()
