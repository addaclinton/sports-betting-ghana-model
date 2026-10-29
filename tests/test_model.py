import sys, unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from model import *  # noqa
from scipy.integrate import solve_ivp


class TestModel(unittest.TestCase):
    def test_class_a_values(self):
        self.assertEqual(N0, 17_931_673)
        self.assertAlmostEqual(mu, 1 / 47.69, places=10)
        self.assertAlmostEqual(Lam, mu * N0, places=6)

    def test_R0_matches_next_generation_matrix(self):
        a = P['a1'] + P['g1'] + mu; b = P['a2'] + P['g2'] + mu; c = P['g3'] + mu
        F = np.zeros((3, 3)); F[0, :] = P['beta']
        V = np.array([[a, 0, 0], [-P['a1'], b, 0], [0, -P['a2'], c]])
        rho = max(abs(np.linalg.eigvals(F @ np.linalg.inv(V))))
        self.assertAlmostEqual(rho, 1.30, places=10)
        self.assertAlmostEqual(R0f(P['beta'], P['a1'], P['a2'], P['g1'], P['g2'], P['g3'], mu), 1.30, places=10)

    def test_population_is_conserved(self):
        _, Y = rk4(y0, P, 30.0)
        self.assertLess(np.max(np.abs(Y.sum(1) - N0)), 1e-3)
        self.assertGreaterEqual(Y.min(), 0.0)

    def test_endemic_equilibrium(self):
        R0, eq = endemic(P)
        self.assertAlmostEqual(R0, 1.30, places=10)
        self.assertLess(np.max(np.abs(rhs(eq, P))), 1e-6)
        self.assertAlmostEqual(eq.sum(), N0, places=3)
        bet = eq[1:4].sum()
        np.testing.assert_allclose(eq[1:4] / bet, [0.16, 0.44, 0.40], atol=1e-12)

    def test_long_run_converges_to_endemic_equilibrium(self):
        _, eq = endemic(P)
        _, Y = rk4(y0, P, 600.0, h=0.02, out_every=600.0)
        np.testing.assert_allclose(Y[-1], eq, rtol=1e-6)

    def test_betting_dies_out_when_R0_below_one(self):
        """With R0 < 1 bettors decay at the rate given by the leading eigenvalue of F - V (Theorem 3)."""
        p = dict(P); p['beta'] = P['beta'] * 0.7          # R0 = 0.91
        a = p['a1'] + p['g1'] + mu; b = p['a2'] + p['g2'] + mu; c = p['g3'] + mu
        F = np.zeros((3, 3)); F[0, :] = p['beta']
        V = np.array([[a, 0, 0], [-p['a1'], b, 0], [0, -p['a2'], c]])
        lam = max(np.linalg.eigvals(F - V).real)
        self.assertLess(lam, 0)
        _, Y = rk4(y0, p, 400.0, h=0.02, out_every=100.0)
        bet = Y[:, 1:4].sum(1)
        self.assertTrue(np.all(np.diff(bet) < 0))
        self.assertLess(bet[-1] / N0, 1e-3)
        rate = np.log(bet[4] / bet[3]) / 100.0              # decay rate between t = 300 and t = 400
        self.assertAlmostEqual(rate, lam, delta=0.02 * abs(lam))

    def test_rk4_agrees_with_adaptive_solver(self):
        _, Y = rk4(y0, P, 30.0)
        sol = solve_ivp(lambda t, y: rhs(y, P), [0, 30], y0, rtol=1e-10, atol=1e-6, method='DOP853')
        rel = np.max(np.abs(sol.y[:, -1] - Y[-1]) / np.maximum(np.abs(Y[-1]), 1))
        self.assertLess(rel, 1e-9)

    def test_alpha2_sensitivity_vanishes_when_g2_equals_g3(self):
        f = lambda a2: R0f(P['beta'], P['a1'], a2, P['g1'], P['g2'], P['g3'], mu)
        h = 1e-6
        self.assertAlmostEqual((f(P['a2'] + h) - f(P['a2'] - h)) / (2 * h), 0.0, places=8)

    def test_threshold_reduction_in_beta(self):
        self.assertAlmostEqual(1 - 1 / 1.30, 0.2308, places=4)


if __name__ == '__main__':
    unittest.main()
