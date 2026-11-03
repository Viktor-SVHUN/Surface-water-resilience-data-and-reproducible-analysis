"""Scientific smoke tests and regression checks against distributed tables."""
from pathlib import Path
import importlib.util
import sys
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('surface_water_resilience', ROOT/'src/surface_water_resilience.py')
ANALYSIS = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ANALYSIS
SPEC.loader.exec_module(ANALYSIS)


class AnalysisTests(unittest.TestCase):
    def test_algorithmic_self_tests(self):
        tests = ANALYSIS.self_tests()
        self.assertEqual(len(tests), 24)
        self.assertTrue(tests.passed.all(), tests.to_string(index=False))

    def test_input_discovery(self):
        self.assertEqual(ANALYSIS.discover_input(), ROOT/'data/raw')

    def test_calendar_gap_and_censoring(self):
        result = ANALYSIS.local_summary([2000, 2001, 2004], [2., 2., 2.], 1., 1.)
        self.assertEqual(result['Lmax_bad_years'], 2)
        self.assertEqual(result['recovery_status'], 'right_censored_no_observed_recovery')

    def test_visibility_monotonic_invariance(self):
        x = np.array([2., 1., 4., 1., 3.])
        self.assertEqual(ANALYSIS.hvg_edges(x), ANALYSIS.hvg_edges(np.exp(x)))

    def test_distributed_ranking(self):
        local = pd.read_csv(ROOT/'results/tables/local_resilience_indicators.csv')
        nodes = pd.read_csv(ROOT/'results/tables/hydrological_nodes.csv')
        expected = pd.read_csv(ROOT/'results/tables/water_body_ranking.csv')
        actual = ANALYSIS.ranking(local, nodes)
        joined = expected.merge(actual, on='VOR', suffixes=('_expected', '_actual'), validate='one_to_one')
        np.testing.assert_allclose(joined.score_expected, joined.score_actual, rtol=1e-12, atol=1e-12)
        self.assertTrue((joined.critical_expected == joined.critical_actual).all())
        self.assertEqual(int(actual.eligible.sum()), 258)
        self.assertEqual(int(actual.critical.sum()), 26)

    def test_clean_source(self):
        source = (ROOT/'src/surface_water_resilience.py').read_text(encoding='utf-8')
        for forbidden in ['D:\\', 'reviewer', 'TODO', 'FIXME', '\\chk{']:
            self.assertNotIn(forbidden, source)


if __name__ == '__main__':
    unittest.main()
