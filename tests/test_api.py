import numpy as np
import pytest
from sklearn.base import clone
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LinearRegression, LogisticRegression

from warrantxai.api import WarrantIndexEvaluator
from warrantxai.dgp import gaussian_linear, gaussian_logistic


def fixture_models():
    dgp = gaussian_linear([2, -1, 0.5], noise=0)
    models = []
    data = []
    for site in range(2):
        replicas = []
        for rep in range(2):
            x, y = dgp.sample(200, site * 10 + rep)
            replicas.append(LinearRegression().fit(x, y))
        models.append(replicas)
        data.append(dgp.sample(1000, 100 + site)[0])
    return models, data


def test_sklearn_clone_end_to_end_and_report():
    models, data = fixture_models()
    estimator = WarrantIndexEvaluator(claim_features=[0, 1], a=0.5)
    cloned = clone(estimator)
    assert cloned.get_params()["a"] == 0.5
    assert cloned.set_params(eta=0.01) is cloned
    assert cloned.fit(models, data) is cloned
    result = cloned.evaluate()
    np.testing.assert_allclose(result["components"][:, 0], 1)
    assert result["WI_claim"] > 0.9
    report = cloned.feature_report()
    assert len(report) == 3
    assert report.point_adequate.all()
    assert (report.certification.str.startswith("not performed")).all()
    result["components"][:] = 0
    assert cloned.components_.max() == 1  # defensive copy


def test_custom_wrong_attributor_and_missing_support():
    models, data = fixture_models()

    def wrong(model, anchors, bg):
        return -(anchors - bg) * model.coef_

    estimator = WarrantIndexEvaluator(attributor=wrong).fit(models, data)
    assert estimator.components_[:, 0].max() < 0.5
    missing = WarrantIndexEvaluator(common_support=False).fit(models, data)
    report = missing.feature_report()
    assert report.WI.isna().all()
    assert report.point_adequate.isna().all()
    single = WarrantIndexEvaluator().fit(models[:1], data[:1])
    assert np.isnan(single.evaluate()["WI"]).all()


def test_binary_logit_adapter():
    dgp = gaussian_logistic([1, -1, 0.5])
    x, y = dgp.sample(500, 77)
    model = LogisticRegression().fit(x, y)
    e = WarrantIndexEvaluator(anchors=x[:100], background_mean=[0, 0, 0]).fit([[model, model], [model, model]], [x, x])
    expected = x[:100] @ model.coef_[0]
    np.testing.assert_allclose(e.reference_[0, 0].sum(1), expected)


def test_not_fitted_and_validation():
    with pytest.raises(NotFittedError):
        WarrantIndexEvaluator().evaluate()
    models, data = fixture_models()
    bad_model = LinearRegression()
    multi = LinearRegression()
    multi.coef_ = np.ones((2, 3))
    calls = [
        lambda: WarrantIndexEvaluator().fit([], []),
        lambda: WarrantIndexEvaluator().fit(models, [data[0], data[1][:, :2]]),
        lambda: WarrantIndexEvaluator().fit([[models[0][0]], models[1]], data),
        lambda: WarrantIndexEvaluator(background_mean=[0, 0]).fit(models, data),
        lambda: WarrantIndexEvaluator().fit([[bad_model, bad_model]] * 2, data),
        lambda: WarrantIndexEvaluator().fit([[multi, multi]] * 2, data),
        lambda: WarrantIndexEvaluator(attributor=lambda *args: np.zeros((1, 3))).fit(models, data),
        lambda: WarrantIndexEvaluator(claim_features=[]).fit(models, data),
    ]
    for call in calls:
        with pytest.raises(ValueError):
            call()


def test_failed_refit_invalidates_previous_fit():
    models, data = fixture_models()
    e = WarrantIndexEvaluator().fit(models, data)
    with pytest.raises(ValueError):
        e.fit([], [])
    with pytest.raises(NotFittedError):
        e.evaluate()


def test_signed_claim_default_and_labeled_legacy_sensitivity():
    models, data = fixture_models()
    for model in models[1]:
        model.coef_[0] *= -1
    claim = WarrantIndexEvaluator().fit(models, data)
    legacy = clone(WarrantIndexEvaluator(generalizability_mode="magnitude_only")).fit(models, data)
    assert claim.evaluate()["WI"][0] < 0.2
    assert legacy.evaluate()["WI"][0] > 0.9
    np.testing.assert_allclose(claim.components_[:, [0, 1, 2, 4]], legacy.components_[:, [0, 1, 2, 4]])
    assert "G_claim" in claim.feature_report()
    assert "G_magnitude_only_sensitivity" in legacy.feature_report()
    assert legacy.evaluate()["generalizability_mode"] == "magnitude_only"
