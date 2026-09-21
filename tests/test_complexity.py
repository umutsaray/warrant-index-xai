import numpy as np
import pytest

from warrantxai.components import importance, rank_reproducibility, stability
from warrantxai.complexity import (
    AttributionCache,
    incomplete_rank,
    incomplete_stability,
    sample_pairs,
    stratified_anchor_subsample,
)


def test_full_pair_budgets_reproduce_exact_formulas():
    rng = np.random.default_rng(33)
    phi = rng.normal(size=(2, 5, 20, 7))
    np.testing.assert_allclose(incomplete_stability(phi, 2, 10, 3), stability(phi, 2))
    np.testing.assert_allclose(
        incomplete_rank(importance(phi), 0.02, 10, 6, 3), rank_reproducibility(importance(phi), 0.02)
    )


def test_pairs_cover_triangle_without_duplicates():
    pairs = sample_pairs(8, 28, 7)
    assert np.all(pairs[:, 0] < pairs[:, 1])
    assert len(np.unique(pairs, axis=0)) == 28
    assert pairs.min() == 0 and pairs.max() == 7


def test_incomplete_site_estimator_monte_carlo_matches_exact():
    # One site removes min-selection bias, so both estimates should be unbiased.
    phi = np.random.default_rng(33).normal(size=(1, 5, 10, 6))
    imp = importance(phi)
    s = np.mean([incomplete_stability(phi, 2, 3, seed) for seed in range(400)], axis=0)
    r = np.mean([incomplete_rank(imp, 0.02, 3, 2, seed) for seed in range(400)], axis=0)
    np.testing.assert_allclose(s, stability(phi, 2), atol=0.025)
    np.testing.assert_allclose(r, rank_reproducibility(imp, 0.02), atol=0.045)


def test_stratified_mean_weights_restore_original_population():
    labels = np.array(["small"] + ["large"] * 99)
    indices, weights = stratified_anchor_subsample(labels, 4, 2)
    assert set(labels[indices]) == {"small", "large"}
    assert len(np.unique(indices)) == 4
    assert weights.sum() == pytest.approx(1)
    values = (labels == "small").astype(float) * 100
    assert np.dot(weights, values[indices]) == pytest.approx(values.mean())


def test_cache_roundtrip_and_semantic_invalidation(tmp_path):
    cache = AttributionCache(tmp_path)
    values = np.arange(12).reshape(2, 2, 3)
    key = cache.key(values, metadata={"background": 0, "version": "test"})
    assert cache.get(key) is None
    cache.put(key, values)
    np.testing.assert_array_equal(cache.get(key), values)
    assert key != cache.key(values, metadata={"background": 1, "version": "test"})
    assert key != cache.key(values + 1, metadata={"background": 0, "version": "test"})
    with pytest.raises(ValueError):
        cache.get("../escape")


@pytest.mark.parametrize(
    "call",
    [
        lambda: sample_pairs(2, 2),
        lambda: sample_pairs(1, 0),
        lambda: incomplete_rank(np.ones((2, 2, 3)), 0, 1, 3),
        lambda: stratified_anchor_subsample([1, 2], 1),
        lambda: stratified_anchor_subsample([], 0),
        lambda: stratified_anchor_subsample([1, np.nan], 2),
    ],
)
def test_invalid_sampling_inputs(call):
    with pytest.raises(ValueError):
        call()
