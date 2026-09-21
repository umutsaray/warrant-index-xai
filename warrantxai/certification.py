"""Section 2.4: conditional validity, intersection-union and BH/BY.

These functions do not manufacture p-values from WI scores. BH requires
independence/PRDS; select BY for arbitrary feature dependence. Neither method
repairs invalid component p-values.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from ._validation import FloatArray, unit_array


@dataclass(frozen=True)
class CertificationResult:
    """Section 2.4 outputs; FDR is not a per-feature truth probability."""

    joint_p: FloatArray
    adjusted_p: FloatArray
    rejected: np.ndarray
    method: str
    alpha: float


def intersection_union(pvalues: ArrayLike) -> FloatArray:
    """Return max_k p_jk for all-required adequacy (Section 2.4 / D.1).

    Parameters
    ----------
    pvalues : array_like of shape (..., required_tests)
        Valid one-sided p-values. Component/site dependence is unrestricted.

    Returns
    -------
    ndarray
        Joint p-values; a union null requires rejecting every component null.
    """
    array = unit_array(pvalues, "pvalues")
    if array.ndim < 1 or not array.shape[-1]:
        raise ValueError("pvalues need a nonempty required-tests axis")
    return array.max(axis=-1)


def fdr_correction(pvalues: ArrayLike, alpha: float = 0.05, method: str = "bh") -> tuple[np.ndarray, FloatArray]:
    """Apply BH (default) or BY, retaining original order (Section 2.4 / D.3).

    Parameters
    ----------
    pvalues : array_like of shape (m,)
        Valid p-values for the preregistered full claim family. NaN rejected.
    alpha : float
        Desired FDR level in (0,1).
    method : {'bh', 'by'}
        BH needs independence/PRDS, BY permits arbitrary dependence.

    Returns
    -------
    rejected, adjusted : ndarray
        Rejection mask and monotone adjusted p-values in input order.
    """
    array = unit_array(pvalues, "pvalues")
    if array.ndim != 1 or not 0 < alpha < 1 or method not in {"bh", "by"}:
        raise ValueError("need p-value vector, 0 < alpha < 1, method 'bh'/'by'")
    m = array.size
    order = np.argsort(array, kind="stable")
    factor = np.sum(1 / np.arange(1, m + 1)) if method == "by" else 1.0
    raw = array[order] * m * factor / np.arange(1, m + 1)
    corrected = np.minimum(1.0, np.minimum.accumulate(raw[::-1])[::-1])
    adjusted = np.empty(m)
    adjusted[order] = corrected
    return adjusted <= alpha, adjusted


def certify(pvalues: ArrayLike, alpha: float = 0.05, method: str = "bh") -> CertificationResult:
    """Perform feature-level IUT then BH/BY (Section 2.4).

    Parameters
    ----------
    pvalues : array_like of shape (features, required_tests)
        Valid p-values; caller owns calibration and family definition.
    alpha : float
        FDR level.
    method : {'bh', 'by'}
        Correction type; not inferred from observed correlations.

    Returns
    -------
    CertificationResult
        Supported all-required claims, not point-estimate adequacy.
    """
    array = unit_array(pvalues, "pvalues")
    if array.ndim != 2:
        raise ValueError("certify requires (features,required_tests)")
    joint = intersection_union(array)
    rejected, adjusted = fdr_correction(joint, alpha, method)
    return CertificationResult(joint, adjusted, rejected, method, alpha)


def independent_block_pvalues(block_scores: ArrayLike, tau: ArrayLike) -> FloatArray:
    """Compute Hoeffding p-values for independent audit blocks (Section 2.4 / D.2).

    assumption: blocks are iid end-to-end studies, not bootstrap replicas of
    one dataset. Target is expected fixed-budget block score, not E[min z].

    Parameters
    ----------
    block_scores : array_like of shape (N, features, criteria)
        Independent bounded study scores.
    tau : array_like
        Fixed scalar or one-dimensional criterion thresholds.

    Returns
    -------
    ndarray of shape (features,criteria)
        Super-uniform p-values for E[score] <= tau.
    """
    array = unit_array(block_scores, "block_scores")
    thresholds = unit_array(tau, "tau")
    if array.ndim != 3 or thresholds.ndim > 1 or (thresholds.ndim == 1 and thresholds.shape != (array.shape[-1],)):
        raise ValueError("expected (N,features,criteria) and scalar/criterion thresholds")
    gap = np.maximum(0.0, array.mean(axis=0) - thresholds)
    return np.exp(-2 * array.shape[0] * gap**2)
