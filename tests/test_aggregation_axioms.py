import numpy as np
import pytest
from hypothesis import given, strategies as st

from warrantxai.aggregation import adequate, claim_index, choquet, geometric_mean, veto_gated_sum, warrant_index

unit = st.floats(min_value=0, max_value=1, allow_nan=False, allow_infinity=False)
vectors = st.lists(unit, min_size=5, max_size=5)


@given(unit)
def test_min_idempotent(c):
    assert warrant_index([c] * 5) == c


@given(vectors, vectors)
def test_min_bounded_monotone(x, increments):
    x = np.array(x)
    y = x + (1 - x) * np.array(increments)
    assert 0 <= warrant_index(x) <= 1
    assert warrant_index(x) <= warrant_index(y)
    assert np.all(warrant_index(x) <= x)


def test_manuscript_counterexample():
    z = [1, 1, 1, 0.1, 1]
    assert warrant_index(z) == 0.1
    assert geometric_mean(z) == pytest.approx(0.6309573444801932)
    assert geometric_mean(z) > 0.6


@given(vectors)
def test_choquet_unanimity_is_min(z):
    nu = np.zeros(32)
    nu[-1] = 1
    assert choquet(z, nu) == pytest.approx(warrant_index(z))


@given(vectors)
def test_additive_choquet_is_mean(z):
    nu = np.array([mask.bit_count() / 5 for mask in range(32)])
    assert choquet(z, nu) == pytest.approx(np.mean(z))


def test_veto_boundary_and_non_idempotency():
    assert veto_gated_sum([0.4] * 5, 0.5) == 0
    assert veto_gated_sum([0.5] * 5, 0.5) == pytest.approx(0.5)
    assert adequate([0.5] * 5, 0.5) == 0  # strict, unlike veto
    assert adequate([0.51] * 5, 0.5) == 1


def test_missing_not_imputed():
    z = [1, 1, 1, np.nan, 1]
    for result in [warrant_index(z), geometric_mean(z), veto_gated_sum(z, 0.5), adequate(z)]:
        assert np.isnan(result)
    assert np.isnan(claim_index([0.9, np.nan], [0, 1]))
    assert claim_index([0.9, np.nan], [0]) == 0.9


@pytest.mark.parametrize("z", [[1] * 4, [2] * 5, [np.inf] * 5, []])
def test_bad_scores(z):
    with pytest.raises(ValueError):
        warrant_index(z)


@pytest.mark.parametrize("weights", [[0.2] * 4, [0, 0, 0, 0, 1], [0.3] * 5])
def test_bad_weights(weights):
    with pytest.raises(ValueError):
        geometric_mean([1] * 5, weights)


def test_invalid_capacity_claim_and_threshold():
    nu = np.ones(32)
    with pytest.raises(ValueError):
        choquet([1] * 5, nu)
    nu[0] = 0
    nu[3] = 0
    with pytest.raises(ValueError):
        choquet([1] * 5, nu)
    for indices in [[], [-1], [3], [0, 0], [0.5]]:
        with pytest.raises(ValueError):
            claim_index([0.5, 0.9], indices)
    for function in [adequate, veto_gated_sum]:
        with pytest.raises(ValueError):
            function([1] * 5, [0.5, 0.5])


def test_geometric_zero_and_vectorization():
    z = np.array([[0, 1, 1, 1, 1], [0.5] * 5])
    np.testing.assert_allclose(geometric_mean(z), [0, 0.5])
    np.testing.assert_allclose(warrant_index(z), [0, 0.5])
