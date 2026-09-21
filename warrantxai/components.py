"""Section 2.2 F/S/U/G_claim/R, with Phase 2.1 signed generalizability.

Attributions have axes (site, training replicate, anchor, feature).
Site zero is the source. No function calibrates tolerances on input data.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import norm

from ._validation import FloatArray, finite_array, tolerance

COMPONENT_NAMES = ("F", "S", "U", "G", "R")


def importance(phi: ArrayLike) -> FloatArray:
    """Compute I_dtj = mean_i |phi_bar_dtj| (Section 2.2).

    Parameters
    ----------
    phi : array_like of shape (sites, T, n, p)
        Signed attributions already averaged over explainer randomness.

    Returns
    -------
    ndarray of shape (sites, T, p)
        Mean absolute importance, not absolute mean attribution.
    """
    return np.abs(finite_array(phi, "phi", 4)).mean(axis=2)


def faithfulness(phi: ArrayLike, reference: ArrayLike, a: ArrayLike) -> FloatArray:
    """Implement F exactly, Section 2.2 F component / Phase 1 B.2.

    Parameters
    ----------
    phi, reference : array_like of shape (sites, T, n, p)
        Candidate phi_bar and independent exact/high-budget psi_tilde.
    a : array_like
        Fixed positive feature tolerance(s).

    Returns
    -------
    ndarray of shape (p,)
        min_site mean_training (1 - mean_anchor clipped squared error).
    """
    x, y = finite_array(phi, "phi", 4), finite_array(reference, "reference", 4)
    if x.shape != y.shape:
        raise ValueError("phi and reference must have identical shapes")
    scale = tolerance(a, x.shape[-1], "a")
    # Clip before squaring to avoid overflow; algebraically identical min(1,r^2).
    loss = np.minimum(1.0, np.abs((x - y) / scale)) ** 2
    return (1 - loss.mean(axis=2)).mean(axis=1).min(axis=0)


def stability(phi: ArrayLike, a: ArrayLike) -> FloatArray:
    """Implement S over all unordered training pairs (Section 2.2 / B.3).

    Parameters
    ----------
    phi : array_like of shape (sites, T, n, p)
        Common-anchor phi_bar; T must be at least two.
    a : array_like
        Fixed positive tolerance(s).

    Returns
    -------
    ndarray of shape (p,)
        Worst-site training stability. Pairs are not independent samples.
    """
    x = finite_array(phi, "phi", 4)
    t = x.shape[1]
    if t < 2:
        raise ValueError("S needs at least two training replicates")
    scale = tolerance(a, x.shape[-1], "a")
    loss = np.zeros((x.shape[0], x.shape[-1]))
    for first, second in combinations(range(t), 2):
        loss += (np.minimum(1.0, np.abs((x[:, first] - x[:, second]) / scale)) ** 2).mean(axis=1)
    return (1 - loss / (t * (t - 1) / 2)).min(axis=0)


def uncertainty(lower: ArrayLike, upper: ArrayLike, u: ArrayLike) -> FloatArray:
    """Implement U from conditional importance CI widths (Section 2.2 / B.4).

    Parameters
    ----------
    lower, upper : array_like of shape (sites, T, p)
        CI limits; upper >= lower. Coverage is a separate obligation.
    u : array_like
        Fixed positive width tolerance(s).

    Returns
    -------
    ndarray of shape (p,)
        High values mean narrow intervals, not calibrated credibility.
    """
    lo, hi = finite_array(lower, "lower", 3), finite_array(upper, "upper", 3)
    if lo.shape != hi.shape or np.any(hi < lo):
        raise ValueError("CI limits must have identical shapes with upper >= lower")
    scale = tolerance(u, lo.shape[-1], "u")
    return (1 - np.minimum(1.0, (hi - lo) / scale)).mean(axis=1).min(axis=0)


def generalizability(site_importance: ArrayLike, b: ArrayLike, common_support: ArrayLike = True) -> FloatArray:
    """Implement magnitude-only G with explicit support mask (Section 2.2 / B.5).

    Parameters
    ----------
    site_importance : array_like of shape (sites, p)
        I_hat_dj, on common anchors/background; site 0 is the source.
    b : array_like
        Fixed positive tolerance(s).
    common_support : bool or array_like of bool
        Scalar, (external_sites,), or (external_sites,p). All preregistered
        external sites must be assessable. No site is silently discarded.

    Returns
    -------
    ndarray of shape (p,)
        Worst-site G; NaN for absent external data or unsupported features.
    """
    imp = finite_array(site_importance, "site_importance", 2)
    if np.any(imp < 0):
        raise ValueError("importance must be nonnegative")
    scale = tolerance(b, imp.shape[-1], "b")
    external, p = imp.shape[0] - 1, imp.shape[-1]
    if not external:
        return np.full(p, np.nan)
    support = np.asarray(common_support)
    if support.dtype.kind != "b":
        raise ValueError("common_support must contain booleans")
    if support.shape == (external,):
        support = support[:, None]
    try:
        support = np.broadcast_to(support, (external, p))
    except ValueError as exc:
        raise ValueError("common_support shape must match external sites/features") from exc
    scores = 1 - np.minimum(1.0, np.abs(imp[1:] - imp[0]) / scale)
    return np.where(support.all(axis=0), scores.min(axis=0), np.nan)


def generalizability_signed(phi_hat: ArrayLike, a: ArrayLike, common_support: ArrayLike = True) -> FloatArray:
    """Implement G_signed exactly (corrected Section 2.2 / Phase 1 B.5).

    Parameters
    ----------
    phi_hat : array_like of shape (sites, n, p)
        Signed site estimates on identical anchors/background. For an ensemble,
        assumption: phi_hat is the training-replicate mean, BEFORE clipping.
    a : array_like
        Fixed positive feature tolerance(s).
    common_support : bool or array_like of bool
        Scalar, (external_sites,), or (external_sites,p); no sites omitted.

    Returns
    -------
    ndarray of shape (p,)
        min_d [1 - mean_i min(1, abs(phi_hat_di-phi_hat_0i)/a)].
        Missing external support or sites gives NaN, not an extremal score.
    """
    x = finite_array(phi_hat, "phi_hat", 3)
    scale = tolerance(a, x.shape[-1], "a")
    # Reuse the unchanged G support validation; identical zero importances
    # give a neutral one where assessable and NaN elsewhere.
    support = generalizability(np.zeros((x.shape[0], x.shape[-1])), 1, common_support)
    if x.shape[0] == 1:
        return support
    scores = 1 - np.minimum(1.0, np.abs(x[1:] - x[0]) / scale).mean(axis=1)
    return np.minimum(support, scores.min(axis=0))


def pairwise_orders(imp: ArrayLike, eta: float) -> NDArray[np.int8]:
    """Form eta-tolerant signed orders o_jl (Section 2.2 / B.6).

    Parameters
    ----------
    imp : array_like of shape (..., p)
        Nonnegative importance.
    eta : float
        Fixed nonnegative tie tolerance; exact boundary is tied.

    Returns
    -------
    ndarray of int8, shape (..., p, p)
        Values in {-1,0,1}; diagonal zero.
    """
    array = finite_array(imp, "importance")
    if array.ndim < 1 or np.any(array < 0) or not np.isfinite(eta) or eta < 0:
        raise ValueError("nonnegative importance and finite eta >= 0 required")
    diff = array[..., :, None] - array[..., None, :]
    return (diff > eta).astype(np.int8) - (diff < -eta).astype(np.int8)


def rank_reproducibility(imp: ArrayLike, eta: float) -> FloatArray:
    """Implement feature-wise R exactly (Section 2.2 / B.6).

    Parameters
    ----------
    imp : array_like of shape (sites, T, p)
        I_dtj; T >= 2 and p >= 2.
    eta : float
        Fixed tie tolerance.

    Returns
    -------
    ndarray of shape (p,)
        Worst-site agreement. Persistent ties score 1, not rank identification.
    """
    array = finite_array(imp, "importance", 3)
    sites, t, p = array.shape
    if t < 2 or p < 2:
        raise ValueError("R needs T >= 2 and p >= 2")
    order = pairwise_orders(array, eta)
    loss = np.zeros((sites, p))
    for first, second in combinations(range(t), 2):
        loss += (np.abs(order[:, first] - order[:, second]) / 2).sum(axis=-1) / (p - 1)
    return (1 - loss / (t * (t - 1) / 2)).min(axis=0)


def importance_ci(phi: ArrayLike, alpha: float = 0.05) -> tuple[FloatArray, FloatArray]:
    """Estimate conditional mean-|phi| normal CIs (Section 2.2 / B.4).

    assumption: independent identically distributed anchors, fixed fitted
    model, exact or pre-averaged phi_bar; CLT approximation, not finite-sample
    coverage. Does not include training uncertainty or unknown explainer bias.

    Parameters
    ----------
    phi : array_like of shape (sites, T, n, p)
        Pre-averaged attributions on independent anchors, n >= 2.
    alpha : float
        Two-sided interval error level.

    Returns
    -------
    lower, upper : ndarray
        Nonnegative lower bounds, shape (sites,T,p).
    """
    array = finite_array(phi, "phi", 4)
    if not 0 < alpha < 1 or array.shape[2] < 2:
        raise ValueError("CI requires 0 < alpha < 1 and n >= 2")
    values = np.abs(array)
    mean = values.mean(axis=2)
    half = norm.ppf(1 - alpha / 2) * values.std(axis=2, ddof=1) / np.sqrt(array.shape[2])
    return np.maximum(0.0, mean - half), mean + half


def component_scores(
    phi: ArrayLike,
    reference: ArrayLike,
    lower: ArrayLike,
    upper: ArrayLike,
    *,
    a: ArrayLike,
    b: ArrayLike,
    u: ArrayLike,
    eta: float,
    common_support: ArrayLike = True,
    generalizability_mode: str = "claim",
) -> FloatArray:
    """Compute all five frozen components (Section 2.2).

    Parameters
    ----------
    phi, reference : array_like
        Matching (sites,T,n,p) candidate/reference tensors.
    lower, upper : array_like
        Matching (sites,T,p) conditional importance CI limits.
    a, b, u : array_like
        Preregistered, positive tolerances.
    eta : float
        Rank tie tolerance.
    common_support : array_like of bool
        External site support mask.
    generalizability_mode : {"claim", "magnitude_only"}
        Default G_claim=min(G,G_signed); magnitude_only is a sensitivity.

    Returns
    -------
    ndarray of shape (p,5)
        Columns F,S,U,selected G,R; default selected G is G_claim.
        No missing-value imputation. The generic G slot name is retained.
    """
    x = finite_array(phi, "phi", 4)
    expected = (x.shape[0], x.shape[1], x.shape[3])
    if np.shape(lower) != expected or np.shape(upper) != expected:
        raise ValueError(f"CI limits must have shape {expected}")
    imp = importance(x)
    g = generalizability(imp.mean(axis=1), b, common_support)
    if generalizability_mode == "claim":
        g = np.minimum(g, generalizability_signed(x.mean(axis=1), a, common_support))
    elif generalizability_mode != "magnitude_only":
        raise ValueError("generalizability_mode must be 'claim' or 'magnitude_only'")
    return np.stack(
        (
            faithfulness(x, reference, a),
            stability(x, a),
            uncertainty(lower, upper, u),
            g,
            rank_reproducibility(imp, eta),
        ),
        axis=-1,
    )
