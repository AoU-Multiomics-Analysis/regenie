"""Check the independent quantitative reference without Spark."""
import unittest
import numpy as np
from scipy.stats import t

class ReferenceTests(unittest.TestCase):
    def test_nonzero_offset_uses_residual_phenotype_units(self):
        from scripts.glow.reference_test import quantitative_reference
        n = 60
        cov = np.linspace(-1, 1, n)
        x = np.asarray([(i * 7 % 19) / 10 for i in range(n)])
        y = .4 * x + .2 * cov + np.sin(np.arange(n) * 1.7)
        design = np.column_stack([np.ones(n), cov])
        q = np.linalg.qr(design)[0]
        offset = .1 * np.cos(np.arange(n))
        offset -= q @ (q.T @ offset)
        residual = y - q @ (q.T @ y)
        scale = np.sqrt(residual @ residual / (n - 2))
        response = residual / scale - offset
        xr = x - q @ (q.T @ x)
        beta = xr @ response / (xr @ xr)
        se = np.sqrt(((response - beta * xr) @ (response - beta * xr)) / (n - 3) / (xr @ xr))
        expected = [beta * scale, se * scale, 2 * t.sf(abs(beta / se), n - 3)]
        actual = quantitative_reference(y, cov, x, offset)
        np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-12)
        self.assertGreater(np.linalg.norm(offset), .1)
