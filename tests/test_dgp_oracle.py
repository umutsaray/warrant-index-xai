import numpy as np
import pytest
from scipy.special import expit

from warrantxai.dgp import GaussianDGP, correlated_sem, exact_shapley, gaussian_linear, gaussian_logistic


@pytest.mark.parametrize("rho", [0, 0.5, 0.9, 1])
def test_oracle_importance_matches_monte_carlo(rho):
    dgp = correlated_sem([2, 0, -1, 0.5], rho)
    x, _ = dgp.sample(150_000, 17)
    phi = dgp.oracle_attributions(x)
    expected = np.abs(dgp.beta) * np.sqrt(2 / np.pi)
    np.testing.assert_allclose(dgp.oracle_importance(), expected)
    error = np.abs(np.abs(phi).mean(0) - expected)
    se = np.abs(phi).std(0, ddof=1) / np.sqrt(len(x))
    assert np.all(error <= 5 * se + 1e-12)


def test_finite_background_enumeration_agrees_with_closed_form():
    dgp = gaussian_linear([2, -1, 0.5], mean=[1, 2, 3], sigma=[1, 2, 0.5])
    x, _ = dgp.sample(5, 1)
    bg, _ = dgp.sample(13, 2)
    exact = exact_shapley(lambda a: 4 + a @ dgp.beta, x, bg)
    np.testing.assert_allclose(exact, dgp.oracle_attributions(x, bg.mean(0)), atol=1e-12)
    np.testing.assert_allclose(exact.sum(1), (x - bg.mean(0)) @ dgp.beta, atol=1e-12)


def test_nonlinear_interaction_exact_shapley():
    x = np.array([[2.0, 3.0]])
    exact = exact_shapley(lambda a: a[:, 0] * a[:, 1], x, np.zeros((1, 2)))
    np.testing.assert_allclose(exact, [[3, 3]])


def test_logit_not_probability_and_bernoulli_sampling():
    dgp = gaussian_logistic([1, -2], intercept=0.3)
    x, y = dgp.sample(20000, 24)
    assert set(np.unique(y)) == {0, 1}
    np.testing.assert_allclose(dgp.oracle_attributions(x).sum(1), x @ dgp.beta)
    assert abs(y.mean() - expit(0.3 + x @ dgp.beta).mean()) < 0.02
    prob_phi = exact_shapley(lambda a: expit(0.3 + a @ dgp.beta), x[:4], np.zeros((1, 2)))
    assert not np.allclose(prob_phi, dgp.oracle_attributions(x[:4]))


def test_shifted_folded_normal_importance_and_degenerate_feature():
    dgp = gaussian_linear([2, 1], mean=[1, 3], sigma=[2, 0])
    x, _ = dgp.sample(150000, 12)
    expected = dgp.oracle_importance([0, 0])
    np.testing.assert_allclose(np.abs(dgp.oracle_attributions(x, [0, 0])).mean(0), expected, rtol=0.006)
    assert expected[1] == 3


def test_seed_reproducibility_and_validation():
    dgp = gaussian_linear([1, 2])
    np.testing.assert_array_equal(dgp.sample(10, 2)[0], dgp.sample(10, 2)[0])
    calls = [
        lambda: correlated_sem([1], 0),
        lambda: correlated_sem([1, 2], 1.1),
        lambda: gaussian_linear([1], sigma=-1),
        lambda: gaussian_linear([1], noise=-1),
        lambda: GaussianDGP(np.ones(2), np.ones(2), np.array([[1, 2], [2, 1]])),
        lambda: GaussianDGP(np.ones(2), np.ones(1), np.eye(2)),
        lambda: dgp.sample(0),
        lambda: dgp.oracle_attributions([[1, 2, 3]]),
        lambda: dgp.oracle_importance([1]),
        lambda: exact_shapley(lambda x: x.sum(1), [[1, 2]], [[0, 0]], max_features=1),
        lambda: exact_shapley(lambda x: np.zeros(1), [[1, 2]], [[0, 0], [1, 1]]),
    ]
    for call in calls:
        with pytest.raises(ValueError):
            call()
