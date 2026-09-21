import numpy as np
import pytest
from scipy.stats import false_discovery_control

from warrantxai.certification import certify, fdr_correction, independent_block_pvalues, intersection_union


def test_hand_iut_bh_by():
    p = [[0.001, 0.002], [0.02, 0.04], [0.9, 0.001]]
    result = certify(p)
    np.testing.assert_allclose(result.joint_p, [0.002, 0.04, 0.9])
    np.testing.assert_allclose(result.adjusted_p, [0.006, 0.06, 0.9])
    np.testing.assert_array_equal(result.rejected, [True, False, False])
    _, by = fdr_correction(result.joint_p, method="by")
    np.testing.assert_allclose(by, [0.011, 0.11, 1])


@pytest.mark.parametrize("method", ["bh", "by"])
def test_correction_matches_independent_scipy_implementation(method):
    p = np.random.default_rng(14).uniform(size=120)
    p[:4] = [0, 1, 0.05, 0.05]
    _, actual = fdr_correction(p, method=method)
    np.testing.assert_allclose(actual, false_discovery_control(p, method=method))


def test_iut_type_i_simulation_under_union_null():
    # 200k genuinely independent synthetic experiments. Least-favourable null:
    # one valid Uniform p and four alternative p=0. Also dependent/all-null cases.
    rng = np.random.default_rng(20260916)
    alpha, n = 0.05, 200_000
    uniform = rng.uniform(size=n)
    one_null = np.zeros((n, 5))
    one_null[:, 0] = uniform
    empirical = np.mean(intersection_union(one_null) <= alpha)
    assert empirical <= alpha  # deterministic seed; not a proof by Monte Carlo
    assert abs(empirical - alpha) < 0.002
    dependent = np.repeat(uniform[:, None], 5, axis=1)
    assert np.mean(intersection_union(dependent) <= alpha) <= alpha
    assert np.mean(intersection_union(rng.uniform(size=(n, 5))) <= alpha) < 0.001


def test_block_p_values_and_type_i_bound():
    rng = np.random.default_rng(19)
    # 10k independent studies with 50 independent Bernoulli blocks per study.
    blocks = rng.binomial(1, 0.5, size=(50, 10_000, 1))
    p = independent_block_pvalues(blocks, 0.5)
    assert np.mean(p <= 0.05) <= 0.05
    np.testing.assert_allclose(independent_block_pvalues(np.zeros((4, 2, 5)), 0.5), 1)
    np.testing.assert_allclose(independent_block_pvalues(np.ones((4, 2, 5)), 0.5), np.exp(-2))


@pytest.mark.parametrize(
    "call",
    [
        lambda: intersection_union(0.5),
        lambda: intersection_union([[np.nan]]),
        lambda: fdr_correction([0.1], method="bad"),
        lambda: fdr_correction([[0.1]]),
        lambda: fdr_correction([0.1], alpha=0),
        lambda: certify([0.1, 0.2]),
        lambda: independent_block_pvalues([[0.5]], 0.5),
        lambda: independent_block_pvalues(np.ones((2, 3, 5)), [0.2, 0.2]),
    ],
)
def test_invalid_inputs(call):
    with pytest.raises(ValueError):
        call()
