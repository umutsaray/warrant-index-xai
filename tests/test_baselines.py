import numpy as np
import pytest

from warrantxai.baselines import ERI_AXES, eri_linear_stress, eri_mean, eri_min, single_component
from warrantxai.scenarios import oracle_invariant, stable_but_wrong, zero_model


def test_eri_aggregates_are_exact():
    assert eri_min([1, 1, 0.1]) == 0.1
    assert eri_mean([1, 1, 0.1]) == pytest.approx(0.7)
    assert np.isnan(eri_min([1, np.nan]))
    assert np.isnan(eri_mean([1, np.nan]))


def test_actual_stress_axes_are_not_wi_components():
    original = eri_linear_stress(oracle_invariant(n=100), seed=33)
    wrong = eri_linear_stress(stable_but_wrong(n=100), seed=33)
    assert tuple(original) == ERI_AXES
    assert all(0 <= v <= 1 for v in original.values())
    # Sign reversal preserves every distance while WI F distinguishes it.
    np.testing.assert_allclose(list(original.values()), list(wrong.values()))
    assert original["ERI-S"] < 1
    np.testing.assert_allclose(list(eri_linear_stress(zero_model(n=100)).values()), 1)


def test_ablations_and_invalid_inputs():
    z = np.array([[0.1, 0.2, 0.3, 0.4, 0.5]])
    assert single_component(z, "F")[0] == 0.1
    assert single_component(z, "S")[0] == 0.2
    assert single_component(z, "G")[0] == 0.4
    for call in [
        lambda: eri_min(0.5),
        lambda: eri_mean(0.5),
        lambda: single_component(z, "U"),
        lambda: eri_linear_stress(oracle_invariant(n=10), probes=1),
    ]:
        with pytest.raises(ValueError):
            call()
