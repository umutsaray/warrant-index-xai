"""Section 2.3: WI and sensitivity aggregators. NaN always propagates."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from ._validation import FloatArray, unit_array


def _scores(z: ArrayLike) -> FloatArray:
    """Require F/S/U/G_claim/R (or explicit magnitude sensitivity), Section 2.3.

    Returns
    -------
    ndarray
        Scores with last axis length five.
    """
    array = unit_array(z, "z", allow_nan=True)
    if array.ndim < 1 or array.shape[-1] != 5:
        raise ValueError("z: last axis must have length 5 (F,S,U,G,R)")
    return array


def _weights(weights: ArrayLike | None) -> FloatArray:
    """Validate strictly positive normalized weights (Section 2.3).

    Returns
    -------
    ndarray
        Five weights; assumption: equal weights when omitted.
    """
    w = np.full(5, 0.2) if weights is None else unit_array(weights, "weights")
    if w.shape != (5,) or np.any(w <= 0) or not np.isclose(w.sum(), 1, atol=1e-12, rtol=0):
        raise ValueError("weights must be five positive numbers summing to 1")
    return w


def warrant_index(z: ArrayLike) -> FloatArray:
    """Compute WI_j = min_k z_jk, Eq. 2.3.1.

    Bounded, monotone, idempotent, continuous, and fully weakest-link.
    NaN means not assessable; it is never omitted from the minimum.

    Parameters
    ----------
    z : array_like of shape (..., 5)
        F, S, U, G_claim, R scores by default. component_scores/API select
        G_claim or explicit magnitude_only sensitivity before aggregation;
        this function does not infer the meaning of a numeric fourth column.

    Returns
    -------
    ndarray
        WI over the last axis.
    """
    return np.min(_scores(z), axis=-1)


def adequate(z: ArrayLike, tau: ArrayLike = 0.8) -> FloatArray:
    """Compute STRICT point-estimate adequacy, Eq. 2.3.2.

    Parameters
    ----------
    z : array_like of shape (..., 5)
        Component scores; this is not statistical certification.
    tau : array_like
        Fixed scalar or five criterion thresholds.

    Returns
    -------
    ndarray
        1 for all z > tau, 0 otherwise, NaN if any dimension unassessed.
    """
    array = _scores(z)
    thresholds = unit_array(tau, "tau")
    if thresholds.ndim > 1 or (thresholds.ndim == 1 and thresholds.shape != (5,)):
        raise ValueError("tau must be scalar or length 5")
    result = np.all(array > thresholds, axis=-1).astype(float)
    return np.where(np.any(np.isnan(array), axis=-1), np.nan, result)


def claim_index(wi: ArrayLike, claim_features: ArrayLike) -> float:
    """Compute min over a preregistered nonempty J*, Eq. 2.3.3.

    Parameters
    ----------
    wi : array_like of shape (p,)
        Feature WI values.
    claim_features : array_like of int
        Unique, preregistered feature indices; never selected here.

    Returns
    -------
    float
        Claim WI, including NaN when a selected feature is unassessed.
    """
    array = unit_array(wi, "wi", allow_nan=True)
    indices = np.asarray(claim_features)
    if array.ndim != 1 or indices.ndim != 1 or not indices.size or indices.dtype.kind not in "iu":
        raise ValueError("expected wi vector and nonempty integer claim_features")
    if np.any(indices < 0) or np.any(indices >= len(array)) or len(np.unique(indices)) != len(indices):
        raise ValueError("claim_features must be unique valid indices")
    return float(np.min(array[indices]))


def geometric_mean(z: ArrayLike, weights: ArrayLike | None = None) -> FloatArray:
    """Compute weighted geometric mean (Section 2.3 sensitivity only).

    Bounded, monotone, idempotent and continuous, but partially compensatory;
    only zero has an absolute veto. Example (1,1,1,.1,1) gives ~.631.

    Parameters
    ----------
    z : array_like of shape (..., 5)
        Component scores.
    weights : array_like of shape (5,), optional
        Strictly positive weights summing to one.

    Returns
    -------
    ndarray
        Geometric mean, without epsilon substitution at zero.
    """
    array, w = _scores(z), _weights(weights)
    with np.errstate(divide="ignore"):
        return np.exp(np.sum(w * np.log(array), axis=-1))


def veto_gated_sum(z: ArrayLike, veto: ArrayLike, weights: ArrayLike | None = None) -> FloatArray:
    """Compute inclusive veto-gated sum (Section 2.3 sensitivity).

    Bounded/monotone; hard veto below thresholds, compensatory above them.
    Not globally idempotent or continuous. This is ELECTRE-inspired, not
    an implementation of ELECTRE III outranking.

    Parameters
    ----------
    z : array_like of shape (..., 5)
        Component scores.
    veto : array_like
        Scalar or length-five inclusive thresholds (z >= v passes).
    weights : array_like, optional
        Positive normalized weights.

    Returns
    -------
    ndarray
        Gated sum; NaNs propagate even when another criterion vetoes.
    """
    array, w = _scores(z), _weights(weights)
    v = unit_array(veto, "veto")
    if v.ndim > 1 or (v.ndim == 1 and v.shape != (5,)):
        raise ValueError("veto must be scalar or length 5")
    result = np.all(array >= v, axis=-1) * np.sum(array * w, axis=-1)
    return np.where(np.any(np.isnan(array), axis=-1), np.nan, result)


def choquet(z: ArrayLike, capacity: ArrayLike) -> FloatArray:
    """Compute the normalized monotone Choquet integral (Section 2.3).

    Bounded, monotone, idempotent and continuous. Not generally veto-based
    or non-compensatory; unanimity capacity reduces exactly to min.

    Parameters
    ----------
    z : array_like of shape (..., 5)
        Component scores.
    capacity : array_like of shape (32,)
        nu(A) indexed by bitmask; bit 0 is F. nu(empty)=0, nu(full)=1.

    Returns
    -------
    ndarray
        Integral for every score vector.
    """
    array = _scores(z)
    nu = unit_array(capacity, "capacity")
    if nu.shape != (32,) or nu[0] != 0 or nu[-1] != 1:
        raise ValueError("capacity requires 32 entries, nu[0]=0, nu[31]=1")
    for mask in range(32):
        for bit in range(5):
            if nu[mask] > nu[mask | (1 << bit)]:
                raise ValueError("capacity must be monotone under set inclusion")
    order = np.argsort(array, axis=-1, kind="stable")
    sorted_z = np.take_along_axis(array, order, axis=-1)
    increments = np.diff(sorted_z, axis=-1, prepend=0.0)
    bits = np.left_shift(1, order)
    masks = np.flip(np.cumsum(np.flip(bits, axis=-1), axis=-1), axis=-1)
    return np.sum(increments * nu[masks], axis=-1)
