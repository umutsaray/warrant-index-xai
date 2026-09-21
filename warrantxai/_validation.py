"""Input contracts shared by the Section 2 estimators."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]


def finite_array(value: ArrayLike, name: str, ndim: int | None = None) -> FloatArray:
    """Validate finite, nonempty numerical inputs (Section 2.2).

    Parameters
    ----------
    value : array_like
        Input data.
    name : str
        Name for validation errors.
    ndim : int, optional
        Required dimensionality.

    Returns
    -------
    ndarray
        Float64 array.
    """
    array = np.asarray(value, dtype=float)
    if array.size == 0 or (ndim is not None and array.ndim != ndim):
        raise ValueError(f"{name}: expected nonempty array with ndim={ndim}")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name}: all values must be finite")
    return array


def tolerance(value: ArrayLike, p: int, name: str) -> FloatArray:
    """Validate positive fixed tolerances, never infer them (Section 2.2).

    Parameters
    ----------
    value : array_like
        Scalar or vector of p positive calibration constants.
    p : int
        Feature count.
    name : str
        Error label.

    Returns
    -------
    ndarray
        Broadcast tolerance vector.
    """
    array = finite_array(value, name)
    if array.ndim > 1 or (array.ndim == 1 and array.shape != (p,)):
        raise ValueError(f"{name}: expected scalar or ({p},)")
    if np.any(array <= 0):
        raise ValueError(f"{name}: tolerances must be strictly positive")
    return np.broadcast_to(array, (p,))


def unit_array(value: ArrayLike, name: str, allow_nan: bool = False) -> FloatArray:
    """Validate [0,1] quantities; optional NaN is unassessed (Section 2.3).

    Parameters
    ----------
    value : array_like
        Scores or probabilities.
    name : str
        Error label.
    allow_nan : bool
        Whether explicit missing scores are permitted.

    Returns
    -------
    ndarray
        Validated scores without imputation or clipping.
    """
    array = np.asarray(value, dtype=float)
    if array.size == 0 or np.any(np.isinf(array)):
        raise ValueError(f"{name}: empty or infinite values")
    if (not allow_nan and np.any(np.isnan(array))) or np.any((array < 0) | (array > 1)):
        raise ValueError(f"{name}: expected values in [0,1]")
    return array
