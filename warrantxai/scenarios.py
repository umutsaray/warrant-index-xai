"""Eleven parameterized synthetic failure scenarios (Section C.2).

All tolerances are external. These fixtures generate controlled model and
explanation outputs, not trained black-box SHAP integrations. Analytic linear
oracles isolate failure modes; selected scenarios genuinely refit OLS models.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from numpy.typing import ArrayLike
from scipy.stats import norm
from sklearn.linear_model import LinearRegression

from ._validation import FloatArray, finite_array
from .components import importance_ci
from .dgp import correlated_sem

DEFAULT_BETA = (2.0, -1.5, 1.0, -0.5, 0.25)


@dataclass
class Scenario:
    """Data-only audit fixture, Section C.2.

    phi/reference have (sites,T,n,p) axes. Coefficients define explanation
    behaviour at arbitrary inputs for independent ERI stress transformations.
    common_support is an explicit structural assertion, not an inferred test.
    """

    name: str
    anchors: FloatArray
    native_anchors: FloatArray
    model_coefficients: FloatArray
    explanation_coefficients: FloatArray
    phi: FloatArray
    reference: FloatArray
    lower: FloatArray
    upper: FloatArray
    common_support: bool | np.ndarray = True
    assumptions: list[str] = field(default_factory=list)
    diagnostics: dict = field(default_factory=dict)


def _base(
    n: int, t: int, sites: int, seed: int, rho: float, beta: ArrayLike
) -> tuple[np.random.Generator, FloatArray, FloatArray, FloatArray]:
    """Generate matched anchors and oracle model coefficients (Section C.2).

    Returns
    -------
    rng, anchors, native_anchors, coefficients : tuple
        Reusable deterministic fixture ingredients.
    """
    if n < 2 or t < 2 or sites < 2:
        raise ValueError("scenarios require n >= 2, T >= 2, sites >= 2")
    rng = np.random.default_rng(seed)
    dgp = correlated_sem(beta, rho)
    x, _ = dgp.sample(n, rng)
    native = np.stack([dgp.sample(n, rng)[0] for _ in range(sites)])
    coefs = np.broadcast_to(dgp.beta, (sites, t, len(dgp.beta))).copy()
    return rng, x, native, coefs


def _assemble(
    name: str,
    x: FloatArray,
    native: FloatArray,
    model: FloatArray,
    candidate: FloatArray | None = None,
    support: bool = True,
    assumptions: list[str] | None = None,
    diagnostics: dict | None = None,
) -> Scenario:
    """Create exact marginal linear reference and conditional CIs (Section C.2).

    Returns
    -------
    Scenario
        Full attribution fixture with no tolerance selection.
    """
    candidate = model.copy() if candidate is None else candidate
    reference = model[:, :, None, :] * x[None, None, :, :]
    phi = candidate[:, :, None, :] * x[None, None, :, :]
    lo, hi = importance_ci(phi)
    return Scenario(
        name,
        x,
        native,
        model,
        candidate,
        phi,
        reference,
        lo,
        hi,
        support,
        ["assumption: marginal linear output, fixed zero background, common anchor Q; normal conditional CIs."]
        + (assumptions or []),
        diagnostics or {},
    )


def oracle_invariant(
    n: int = 1000, t: int = 6, sites: int = 3, seed: int = 0, rho: float = 0, beta: ArrayLike = DEFAULT_BETA
) -> Scenario:
    """Generate oracle-invariant scenario (Section C.2).

    Parameters
    ----------
    n, t, sites : int
        Anchors, model replicas and total domains (source included).
    seed, rho, beta : int, float, array_like
        Reproducibility, feature correlation and known mechanism.

    Returns
    -------
    Scenario
        Exact attributions with identical oracle mechanisms at every site.
    """
    _, x, native, coefs = _base(n, t, sites, seed, rho, beta)
    return _assemble(
        "oracle_invariant",
        x,
        native,
        coefs,
        assumptions=[
            "assumption: independent populations share a stipulated oracle mechanism; no training estimation noise."
        ],
    )


def stable_but_wrong(
    n: int = 1000,
    t: int = 6,
    sites: int = 3,
    seed: int = 0,
    rho: float = 0,
    beta: ArrayLike = DEFAULT_BETA,
    severity: float = 1,
) -> Scenario:
    """Generate a stable, sign-reversed wrong explainer (Section C.2).

    Parameters
    ----------
    n, t, sites, seed, rho, beta : see oracle_invariant
        Shared DGP parameters.
    severity : float in [0,1]
        Interpolates from correct to sign-reversed; model is unchanged.

    Returns
    -------
    Scenario
        S/R remain high despite low F at full severity.
    """
    if not 0 <= severity <= 1:
        raise ValueError("severity must be in [0,1]")
    _, x, native, coefs = _base(n, t, sites, seed, rho, beta)
    return _assemble("stable_but_wrong", x, native, coefs, coefs * (1 - 2 * severity))


def noisy(
    n: int = 1000,
    t: int = 6,
    sites: int = 3,
    seed: int = 0,
    rho: float = 0,
    beta: ArrayLike = DEFAULT_BETA,
    noise_sd: float = 2,
    mc_draws: int = 4,
) -> Scenario:
    """Generate controlled finite-MC explainer uncertainty (Section C.2).

    Parameters
    ----------
    n, t, sites, seed, rho, beta : see oracle_invariant
        Shared DGP parameters.
    noise_sd : float
        Known coefficient-estimation standard deviation per MC draw.
    mc_draws : int
        Number of independent draws used to estimate phi_bar.

    Returns
    -------
    Scenario
        Finite-MC mean attributions and uncertainty-aware importance CIs.
    """
    if not np.isfinite(noise_sd) or noise_sd < 0 or mc_draws < 1:
        raise ValueError("noise_sd >= 0 and mc_draws >= 1 required")
    rng, x, native, model = _base(n, t, sites, seed, rho, beta)
    errors = rng.normal(0, noise_sd, size=(mc_draws, len(beta)))
    candidate = model + errors.mean(0)
    result = _assemble(
        "noisy",
        x,
        native,
        model,
        candidate,
        assumptions=[
            "assumption: common random numbers across training/site replicas isolate MC noise from training variability.",
            "assumption: known Gaussian coefficient MC error; Bonferroni-split coefficient/anchor normal intervals form magnitude-product CIs; anchor part is asymptotic.",
            "assumption: finite MC estimate of phi_bar, not its infinite-sample expectation.",
        ],
    )
    critical = norm.ppf(1 - 0.05 / 4)
    coefficient_radius = critical * noise_sd / np.sqrt(mc_draws)
    mean_x = np.abs(x).mean(0)
    anchor_radius = critical * np.abs(x).std(0, ddof=1) / np.sqrt(n)
    result.lower = np.maximum(0, np.abs(candidate) - coefficient_radius) * np.maximum(0, mean_x - anchor_radius)
    result.upper = (np.abs(candidate) + coefficient_radius) * (mean_x + anchor_radius)
    result.diagnostics = {"mc_draws": mc_draws, "coefficient_mc_se": noise_sd / np.sqrt(mc_draws)}
    return result


def training_instability(
    n: int = 1000,
    t: int = 6,
    sites: int = 3,
    seed: int = 0,
    rho: float = 0.9,
    beta: ArrayLike = DEFAULT_BETA,
    train_n: int = 30,
    noise_sd: float = 5,
) -> Scenario:
    """Refit noisy small-sample OLS independently (Section C.2).

    Parameters
    ----------
    n, t, sites, seed, rho, beta : see oracle_invariant
        Shared anchor and DGP parameters.
    train_n : int
        Independent training observations per replica.
    noise_sd : float
        Outcome residual standard deviation.

    Returns
    -------
    Scenario
        Faithful explanations of unstable fitted models.
    """
    rng, x, native, coefs = _base(n, t, sites, seed, rho, beta)
    dgp = correlated_sem(beta, rho, noise=noise_sd)
    for d in range(sites):
        for rep in range(t):
            train, y = dgp.sample(train_n, rng)
            coefs[d, rep] = LinearRegression().fit(train, y).coef_
    return _assemble(
        "training_instability",
        x,
        native,
        coefs,
        assumptions=[
            "assumption: independent OLS training samples (not bootstrap replicas) with finite-sample model variation."
        ],
    )


def covariate_shift(
    n: int = 1000,
    t: int = 6,
    sites: int = 3,
    seed: int = 0,
    rho: float = 0,
    beta: ArrayLike = DEFAULT_BETA,
    shift: float = 2,
) -> Scenario:
    """Change population means/scales, preserve the mechanism (Section C.2).

    Parameters
    ----------
    n, t, sites, seed, rho, beta : see oracle_invariant
        Shared parameters and common-Q anchors.
    shift : float
        External site mean shift; finite Gaussian support is shared.

    Returns
    -------
    Scenario
        Standardized G=1 while native-population importance changes.
    """
    if not np.isfinite(shift):
        raise ValueError("shift must be finite")
    rng, x, native, coefs = _base(n, t, sites, seed, rho, beta)
    for d in range(1, sites):
        native[d] = native[d] * (1 + 0.3 * d) + shift * d
    return _assemble(
        "covariate_shift",
        x,
        native,
        coefs,
        assumptions=[
            "assumption: all site models are oracle beta; comparison uses shared Q, not native shifted cohorts."
        ],
        diagnostics={"native_mean_absolute_attribution": np.abs(native * np.asarray(beta)).mean(1).tolist()},
    )


def mechanism_shift(
    n: int = 1000,
    t: int = 6,
    sites: int = 3,
    seed: int = 0,
    rho: float = 0,
    beta: ArrayLike = DEFAULT_BETA,
    feature: int = 0,
    mode: str = "vanish",
) -> Scenario:
    """Remove or reverse a site's true mechanism (Section C.2).

    Parameters
    ----------
    n, t, sites, seed, rho, beta : see oracle_invariant
        Shared parameters.
    feature : int
        Feature whose external coefficient changes.
    mode : {'vanish','sign_flip'}
        sign_flip is a declared blind-spot diagnostic for magnitude-only G.

    Returns
    -------
    Scenario
        Each explanation is faithful to its own site's model.
    """
    if mode not in {"vanish", "sign_flip"} or not 0 <= feature < len(beta):
        raise ValueError("invalid mechanism mode/feature")
    _, x, native, coefs = _base(n, t, sites, seed, rho, beta)
    coefs[1:, :, feature] *= 0 if mode == "vanish" else -1
    return _assemble("mechanism_shift", x, native, coefs, diagnostics={"mode": mode, "changed_feature": feature})


def spurious_shortcut(
    n: int = 1000,
    t: int = 6,
    sites: int = 3,
    seed: int = 0,
    rho: float = 0,
    beta: ArrayLike = DEFAULT_BETA,
    train_n: int = 500,
    shortcut_noise: float = 0.1,
) -> Scenario:
    """Fit source shortcut Z=Y+noise that disappears externally (Section C.2).

    Parameters
    ----------
    n, t, sites, seed, rho, beta : see oracle_invariant
        Last feature is repurposed as a noncausal outcome proxy.
    train_n : int
        Training sample count per model.
    shortcut_noise : float
        Source proxy noise standard deviation.

    Returns
    -------
    Scenario
        Stable source shortcut, faithful explanations and external G failure.
    """
    if train_n < 2 or not np.isfinite(shortcut_noise) or shortcut_noise <= 0:
        raise ValueError("train_n >= 2 and shortcut_noise > 0 required")
    rng, x, native, coefs = _base(n, t, sites, seed, rho, beta)
    p = len(beta)
    causal = finite_array(beta, "beta", 1).copy()
    causal[-1] = 0
    for d in range(sites):
        for rep in range(t):
            train = rng.normal(size=(train_n, p))
            y = train @ causal + rng.normal(size=train_n)
            train[:, -1] = y + rng.normal(0, shortcut_noise, train_n) if d == 0 else rng.normal(size=train_n)
            coefs[d, rep] = LinearRegression().fit(train, y).coef_
        sample = rng.normal(size=(n, p))
        y = sample @ causal + rng.normal(size=n)
        sample[:, -1] = y + rng.normal(0, shortcut_noise, n) if d == 0 else rng.normal(size=n)
        native[d] = sample
    return _assemble(
        "spurious_shortcut",
        x,
        native,
        coefs,
        assumptions=[
            "assumption: final feature is an outcome proxy, not a causal parent; external populations remove this proxy relationship.",
            "assumption: common Q is a fixed standard-normal audit distribution; positive noise gives common mathematical support, not guaranteed practical overlap.",
        ],
        diagnostics={"shortcut_feature": p - 1},
    )


def near_tied_ranks(
    n: int = 1000,
    t: int = 6,
    sites: int = 3,
    seed: int = 0,
    rho: float = 0,
    beta: ArrayLike = DEFAULT_BETA,
    coefficient_sd: float = 0.08,
) -> Scenario:
    """Perturb nearly equal coefficients across replicas (Section C.2).

    Parameters
    ----------
    n, t, sites, seed, rho, beta : see oracle_invariant
        Shared DGP parameters; first two coefficients have magnitude one.
    coefficient_sd : float
        Known training-coefficient perturbation scale.

    Returns
    -------
    Scenario
        High faithfulness, ambiguous pairwise ranking; WI need not reject.
    """
    if not np.isfinite(coefficient_sd) or coefficient_sd < 0:
        raise ValueError("coefficient_sd must be nonnegative")
    rng, x, native, coefs = _base(n, t, sites, seed, rho, beta)
    coefs[:, :, :2] = 1 + rng.normal(0, coefficient_sd, size=(sites, t, 2))
    return _assemble(
        "near_tied_ranks",
        x,
        native,
        coefs,
        assumptions=[
            "assumption: controlled coefficient-perturbation model ensemble; not an OLS uncertainty guarantee."
        ],
    )


def redundant_duplicates(
    n: int = 1000,
    t: int = 6,
    sites: int = 3,
    seed: int = 0,
    rho: float = 0,
    beta: ArrayLike = DEFAULT_BETA,
    group_coefficient: float = 2,
) -> Scenario:
    """Redistribute a fixed coefficient sum across exact duplicates (Section C.2).

    Parameters
    ----------
    n, t, sites, seed, rho, beta : see oracle_invariant
        Shared parameters; features 0 and 1 become exact duplicates.
    group_coefficient : float
        Fixed predictive group effect split differently across models.

    Returns
    -------
    Scenario
        Same on-support predictions; unstable individual feature attribution.
    """
    if not np.isfinite(group_coefficient):
        raise ValueError("group_coefficient must be finite")
    rng, x, native, coefs = _base(n, t, sites, seed, rho, beta)
    x[:, 1] = x[:, 0]
    native[:, :, 1] = native[:, :, 0]
    splits = rng.uniform(size=(sites, t))
    coefs[:, :, 0] = group_coefficient * splits
    coefs[:, :, 1] = group_coefficient * (1 - splits)
    result = _assemble(
        "redundant_duplicates",
        x,
        native,
        coefs,
        assumptions=[
            "assumption: nonidentifiable coefficient splits under exact collinearity; marginal SHAP probes off-manifold coalitions.",
            "assumption: group-level invariance is diagnostic only; WI remains feature-level as specified.",
        ],
    )
    result.diagnostics = {
        "max_group_attribution_difference": float(np.max(np.ptp(result.phi[..., :2].sum(-1), axis=1)))
    }
    return result


def ood(
    n: int = 1000,
    t: int = 6,
    sites: int = 3,
    seed: int = 0,
    rho: float = 0,
    beta: ArrayLike = DEFAULT_BETA,
    gap: float = 2,
) -> Scenario:
    """Create genuinely disjoint bounded supports (Section C.2).

    Parameters
    ----------
    n, t, sites, seed, rho, beta : see oracle_invariant
        rho is unused after uniform-support replacement (explicit assumption).
    gap : float
        Positive gap separating source [-1,1] and external support.

    Returns
    -------
    Scenario
        G/WI are not assessable, even though the linear oracle can extrapolate.
    """
    if not np.isfinite(gap) or gap <= 0:
        raise ValueError("gap must be positive")
    rng, x, native, coefs = _base(n, t, sites, seed, rho, beta)
    x = rng.uniform(-1, 1, x.shape)
    native[0] = rng.uniform(-1, 1, x.shape)
    native[1:] = rng.uniform(1 + gap, 3 + gap, native[1:].shape)
    return _assemble(
        "ood",
        x,
        native,
        coefs,
        support=False,
        assumptions=[
            "assumption: bounded uniform supports replace Gaussian data to ensure structural non-overlap; rho is unused."
        ],
    )


def zero_model(
    n: int = 1000, t: int = 6, sites: int = 3, seed: int = 0, rho: float = 0, beta: ArrayLike = DEFAULT_BETA
) -> Scenario:
    """Generate constant model with correct zero explanation (Section C.2).

    Parameters
    ----------
    n, t, sites, seed, rho, beta : see oracle_invariant
        beta sets dimensionality only; all model coefficients are zero.

    Returns
    -------
    Scenario
        WI=1 is expected; it cannot imply an informative or unique ranking.
    """
    _, x, native, coefs = _base(n, t, sites, seed, rho, beta)
    coefs[:] = 0
    return _assemble("zero_model", x, native, coefs, diagnostics={"informative": False})


SCENARIOS: dict[str, Callable[..., Scenario]] = {
    function.__name__: function
    for function in [
        oracle_invariant,
        stable_but_wrong,
        noisy,
        training_instability,
        covariate_shift,
        mechanism_shift,
        spurious_shortcut,
        near_tied_ranks,
        redundant_duplicates,
        ood,
        zero_model,
    ]
}
