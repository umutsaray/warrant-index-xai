import numpy as np
import pytest

from warrantxai.aggregation import warrant_index
from warrantxai.components import component_scores
from warrantxai.scenarios import (
    SCENARIOS,
    mechanism_shift,
    noisy,
    redundant_duplicates,
    spurious_shortcut,
    training_instability,
)


def scores(scenario):
    return component_scores(
        scenario.phi,
        scenario.reference,
        scenario.lower,
        scenario.upper,
        a=0.5,
        b=0.5,
        u=3,
        eta=0.02,
        common_support=scenario.common_support,
    )


@pytest.mark.parametrize("name", list(SCENARIOS))
def test_every_scenario_is_reproducible_and_well_formed(name):
    s = SCENARIOS[name](n=150, t=4, sites=2, seed=15)
    repeated = SCENARIOS[name](n=150, t=4, sites=2, seed=15)
    np.testing.assert_array_equal(s.phi, repeated.phi)
    assert s.phi.shape == (2, 4, 150, 5)
    assert s.lower.shape == (2, 4, 5)
    assert (s.lower <= s.upper).all()
    assert s.assumptions
    z = scores(s)
    finite = z[np.isfinite(z)]
    assert np.all((finite >= 0) & (finite <= 1))
    if name == "ood":
        assert np.isnan(z[:, 3]).all()
        assert np.isnan(warrant_index(z)).all()


def test_stable_wrong_detected_by_f_not_s():
    z = scores(SCENARIOS["stable_but_wrong"](seed=9))
    assert np.max(z[:, 0]) < 0.6
    np.testing.assert_allclose(z[:, [1, 3, 4]], 1)


def test_mechanism_vanish_and_sign_flip_claim_regression():
    vanished = scores(mechanism_shift())
    assert vanished[0, 3] == 0
    flipped = scores(mechanism_shift(mode="sign_flip"))
    np.testing.assert_allclose(flipped[:, [0, 1, 4]], 1)
    # Phase 2.1: default G_claim now detects signed changes despite unchanged magnitudes.
    assert flipped[0, 3] < 0.2
    assert warrant_index(flipped)[0] < 0.8


def test_noisy_more_mc_reduces_uncertainty():
    low = scores(noisy(mc_draws=4, seed=21))
    high = scores(noisy(mc_draws=4000, seed=21))
    assert high[:, 2].mean() > low[:, 2].mean() + 0.4
    np.testing.assert_allclose(low[:, 1], 1)  # common random numbers isolate MC source


def test_training_and_shortcut_are_actual_fitted_models():
    unstable = training_instability(seed=17)
    assert np.std(unstable.model_coefficients[:, :, 0]) > 0.1
    z = scores(unstable)
    np.testing.assert_allclose(z[:, 0], 1)
    assert z[:, 1].min() < 0.7
    shortcut = spurious_shortcut(seed=17)
    assert shortcut.model_coefficients[0, :, -1].mean() > 0.8
    assert abs(shortcut.model_coefficients[1, :, -1].mean()) < 0.2
    assert scores(shortcut)[-1, 3] < 0.2


def test_duplicates_group_invariant_and_zero_model_correct():
    duplicates = redundant_duplicates()
    assert duplicates.diagnostics["max_group_attribution_difference"] < 1e-12
    np.testing.assert_allclose(scores(SCENARIOS["zero_model"]()), 1)


@pytest.mark.parametrize(
    "call",
    [
        lambda: SCENARIOS["oracle_invariant"](n=1),
        lambda: SCENARIOS["stable_but_wrong"](severity=2),
        lambda: noisy(mc_draws=0),
        lambda: SCENARIOS["covariate_shift"](shift=np.inf),
        lambda: mechanism_shift(mode="bad"),
        lambda: spurious_shortcut(shortcut_noise=0),
        lambda: SCENARIOS["near_tied_ranks"](coefficient_sd=-1),
        lambda: redundant_duplicates(group_coefficient=np.inf),
        lambda: SCENARIOS["ood"](gap=0),
    ],
)
def test_scenario_invalid_parameters(call):
    with pytest.raises(ValueError):
        call()
