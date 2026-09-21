import numpy as np
import pytest

from warrantxai.components import (
    component_scores,
    faithfulness,
    generalizability,
    generalizability_signed,
    importance,
    importance_ci,
    pairwise_orders,
    rank_reproducibility,
    stability,
    uncertainty,
)


def test_hand_calculation_three_features_two_sites_two_replicates():
    # Source I_t = (1,2,0),(2,1,0); external repeats (1,2,0).
    phi = np.array([[[[1, 2, 0]], [[2, 1, 0]]], [[[1, 2, 0]], [[1, 2, 0]]]], dtype=float)
    reference = phi.copy()
    reference[0, 0, 0, 0] += 1
    np.testing.assert_allclose(faithfulness(phi, reference, 2), [0.875, 1, 1])
    np.testing.assert_allclose(stability(phi, 2), [0.75, 0.75, 1])
    lo = np.zeros((2, 2, 3))
    hi = np.ones((2, 2, 3)) * 0.2
    np.testing.assert_allclose(uncertainty(lo, hi, 1), [0.8] * 3)
    np.testing.assert_allclose(generalizability(importance(phi).mean(1), 1), [0.5, 0.5, 1])
    np.testing.assert_allclose(rank_reproducibility(importance(phi), 0), [0.5, 0.5, 1])
    z = component_scores(phi, reference, lo, hi, a=2, b=1, u=1, eta=0)
    np.testing.assert_allclose(z[:, 0], [0.875, 1, 1])
    assert z.shape == (3, 5)


def test_squares_clipped_before_averaging_and_worst_site():
    phi = np.zeros((2, 2, 2, 2))
    ref = phi.copy()
    ref[1, :, 0, 0] = 20
    np.testing.assert_allclose(faithfulness(phi, ref, 1), [0.5, 1])
    phi[1, 1, 0, 0] = 20
    np.testing.assert_allclose(stability(phi, 1), [0.5, 1])


def test_ties_boundary_and_persistent_zero():
    orders = pairwise_orders([0, 0.5, 1], 0.5)
    assert orders[0, 1] == 0
    assert orders[0, 2] == -1
    np.testing.assert_allclose(rank_reproducibility(np.zeros((2, 2, 3)), 0), 1)
    # tie-to-strict costs 1/2, p=2.
    np.testing.assert_allclose(rank_reproducibility([[[1, 1], [2, 1]]], 0), [0.5, 0.5])


def test_support_and_no_external_site():
    assert np.isnan(generalizability([[1, 2]], 1)).all()
    assert np.isnan(generalizability([[1, 2], [1, 2]], 1, False)).all()
    np.testing.assert_allclose(generalizability([[1, 2], [1, 2]], 1, [[True, False]]), [1, np.nan])
    assert np.isnan(generalizability([[1, 2], [1, 2], [1, 2]], 1, [True, False])).all()
    with pytest.raises(ValueError):
        generalizability([[1, 2], [1, 2]], 1, [1])
    with pytest.raises(ValueError):
        generalizability([[1, 2], [1, 2]], 1, [True] * 3)


def test_ci_and_abs_mean():
    phi = np.tile(np.array([-1, 1.0])[None, None, :, None], (2, 2, 1, 3))
    np.testing.assert_allclose(importance(phi), 1)
    lo, hi = importance_ci(phi)
    np.testing.assert_allclose(lo, 1)
    np.testing.assert_allclose(hi, 1)
    np.testing.assert_allclose(uncertainty(lo, hi, 1), 1)


@pytest.mark.parametrize("a", [0, -1, np.nan, [1, 2, 3], [[1, 1]]])
def test_invalid_tolerances(a):
    x = np.zeros((2, 2, 2, 2))
    with pytest.raises(ValueError):
        faithfulness(x, x, a)


def test_invalid_shapes_and_values():
    x = np.zeros((2, 2, 2, 2))
    for call in [
        lambda: faithfulness(x, x[:1], 1),
        lambda: stability(x[:, :1], 1),
        lambda: rank_reproducibility(np.zeros((2, 1, 2)), 0),
        lambda: rank_reproducibility(np.zeros((2, 2, 1)), 0),
        lambda: pairwise_orders([1, 2], -1),
        lambda: pairwise_orders([-1, 2], 0),
        lambda: uncertainty(np.ones((2, 2, 2)), np.zeros((2, 2, 2)), 1),
        lambda: importance_ci(x, alpha=1),
        lambda: importance_ci(x[:, :, :1]),
        lambda: generalizability([[-1, 2], [1, 2]], 1),
        lambda: component_scores(x, x, [0], [1], a=1, b=1, u=1, eta=0),
        lambda: importance(np.full_like(x, np.nan)),
        lambda: importance([]),
    ]:
        with pytest.raises(ValueError):
            call()


def test_signed_hand_calculation_clipping_then_anchor_mean_then_site_min():
    # a=(2,1): site1 losses=(.75,.5), site2=(.625,.25).
    x = np.array([[[0, 0], [0, 0]], [[1, 0], [4, 1]], [[0.5, 0.5], [2, 0]]])
    np.testing.assert_allclose(generalizability_signed(x, [2, 1]), [0.25, 0.5])
    # First-feature site1 loss is (.5+1)/2=.75 (not clipped mean difference).
    assert np.isnan(generalizability_signed(x[:1], 1)).all()
    assert np.isnan(generalizability_signed(x, 1, False)).all()
    np.testing.assert_allclose(generalizability_signed(x, [2, 1], [[True, False]] * 2), [0.25, np.nan])
    assert np.isnan(generalizability_signed(x, 1, [True, False])).all()
    for support in ([1, 1], [True] * 4):
        with pytest.raises(ValueError):
            generalizability_signed(x, 1, support)
    with pytest.raises(ValueError):
        generalizability_signed(x, 0)


def test_claim_replica_mean_before_clipping_and_magnitude_compatibility():
    phi = np.zeros((2, 2, 2, 2))
    phi[1, 0] = 2
    phi[1, 1] = -2
    ci = np.zeros((2, 2, 2))
    # Signed replica means cancel: no invented pair-wise signed penalty.
    np.testing.assert_allclose(generalizability_signed(phi.mean(1), 1), 1)
    for mode in ("claim", "magnitude_only"):
        z = component_scores(phi, phi, ci, ci, a=1, b=4, u=1, eta=0, generalizability_mode=mode)
        np.testing.assert_allclose(z[:, 3], 0.5)
    with pytest.raises(ValueError, match="generalizability_mode"):
        component_scores(phi, phi, ci, ci, a=1, b=1, u=1, eta=0, generalizability_mode="bad")
