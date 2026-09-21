"""Synthetic Gaussian DGPs and marginal Shapley oracles (Section C.1).

No SHAP/LIME dependency. Correlated features do not change marginal linear
Shapley values; conditional or causal Shapley would be a different estimand.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import factorial
from typing import Callable

import numpy as np
from numpy.typing import ArrayLike
from scipy.special import expit
from scipy.stats import norm

from ._validation import FloatArray, finite_array


@dataclass(frozen=True)
class GaussianDGP:
    """Gaussian linear/logistic DGP (Section C.1).

    Parameters
    ----------
    beta, mean : ndarray of shape (p,)
        Mechanism coefficients and feature mean.
    covariance : ndarray of shape (p,p)
        Positive semidefinite covariance; singular duplicates are allowed.
    intercept, noise : float
        Linear intercept and nonnegative residual standard deviation.
    logistic : bool
        If true, Y is Bernoulli and oracle output is the logit, not probability.
    """

    beta: FloatArray
    mean: FloatArray
    covariance: FloatArray
    intercept: float = 0.0
    noise: float = 1.0
    logistic: bool = False

    def __post_init__(self) -> None:
        """Validate DGP parameters (Section C.1).

        Returns
        -------
        None
        """
        beta = finite_array(self.beta, "beta", 1).copy()
        mean = finite_array(self.mean, "mean", 1).copy()
        cov = finite_array(self.covariance, "covariance", 2).copy()
        if mean.shape != beta.shape or cov.shape != (len(beta), len(beta)):
            raise ValueError("DGP parameter shapes do not match")
        if not np.allclose(cov, cov.T, atol=1e-12, rtol=0) or np.linalg.eigvalsh(cov).min() < -1e-10:
            raise ValueError("covariance must be symmetric positive semidefinite")
        if not np.isfinite(self.intercept) or not np.isfinite(self.noise) or self.noise < 0:
            raise ValueError("intercept must be finite and noise nonnegative")
        for name, value in [("beta", beta), ("mean", mean), ("covariance", cov)]:
            value.setflags(write=False)
            object.__setattr__(self, name, value)

    def sample(self, n: int, seed: int | np.random.Generator = 0) -> tuple[FloatArray, FloatArray]:
        """Sample features and labels from the known mechanism (Section C.1).

        Parameters
        ----------
        n : int
            Positive observation count.
        seed : int or Generator
            Explicit randomness source.

        Returns
        -------
        x, y : ndarray
            Features (n,p) and linear responses or binary labels (n,).
        """
        if not isinstance(n, (int, np.integer)) or n <= 0:
            raise ValueError("n must be a positive integer")
        rng = np.random.default_rng(seed)
        x = rng.multivariate_normal(self.mean, self.covariance, size=n)
        signal = self.intercept + x @ self.beta
        y = rng.binomial(1, expit(signal)) if self.logistic else signal + rng.normal(0, self.noise, n)
        return x, np.asarray(y, dtype=float)

    def oracle_attributions(self, x: ArrayLike, background_mean: ArrayLike | None = None) -> FloatArray:
        """Return beta_j (x_j-mu_bg,j), marginal game (Section C.1).

        assumption: output is linear signal or logistic logit. This is exact
        for marginal/interventional SHAP even under feature correlation.

        Parameters
        ----------
        x : array_like of shape (n,p)
            Anchor observations.
        background_mean : array_like of shape (p,), optional
            Fixed background mean; defaults to the DGP feature mean.

        Returns
        -------
        ndarray of shape (n,p)
            Exact signed marginal attributions.
        """
        anchors = finite_array(x, "x", 2)
        bg = self.mean if background_mean is None else finite_array(background_mean, "background_mean", 1)
        if anchors.shape[1] != len(self.beta) or bg.shape != self.beta.shape:
            raise ValueError("anchor/background feature counts must match beta")
        return (anchors - bg) * self.beta

    def oracle_importance(self, background_mean: ArrayLike | None = None) -> FloatArray:
        """Return exact mean absolute Gaussian attribution (Section C.1).

        Matching means give |beta_j| sigma_j sqrt(2/pi). Otherwise the exact
        folded-normal mean is used, avoiding hidden re-centering under shift.

        Parameters
        ----------
        background_mean : array_like, optional
            Fixed reference mean, default self.mean.

        Returns
        -------
        ndarray of shape (p,)
            Exact marginal importance under this DGP's native population.
        """
        bg = self.mean if background_mean is None else finite_array(background_mean, "background_mean", 1)
        if bg.shape != self.beta.shape:
            raise ValueError("background_mean must match beta")
        sigma = np.sqrt(np.maximum(0, np.diag(self.covariance)))
        offset = np.abs(self.mean - bg)
        ratio = np.divide(offset, sigma, out=np.zeros_like(offset), where=sigma > 0)
        folded = np.where(
            sigma > 0, sigma * np.sqrt(2 / np.pi) * np.exp(-0.5 * ratio**2) + offset * (2 * norm.cdf(ratio) - 1), offset
        )
        return np.abs(self.beta) * folded


def gaussian_linear(
    beta: ArrayLike,
    *,
    mean: ArrayLike | None = None,
    sigma: ArrayLike | float = 1.0,
    noise: float = 1.0,
    intercept: float = 0.0,
) -> GaussianDGP:
    """Create the independent Gaussian linear DGP (Section C.1).

    Parameters
    ----------
    beta : array_like of shape (p,)
        True coefficients.
    mean, sigma : array_like, optional
        Feature means and nonnegative marginal standard deviations.
    noise, intercept : float
        Residual standard deviation and signal intercept.

    Returns
    -------
    GaussianDGP
        Validated independent-feature DGP.
    """
    coefficients = finite_array(beta, "beta", 1)
    mu = np.zeros_like(coefficients) if mean is None else finite_array(mean, "mean", 1)
    scale = np.broadcast_to(finite_array(sigma, "sigma"), coefficients.shape)
    if np.any(scale < 0):
        raise ValueError("sigma must be nonnegative")
    return GaussianDGP(coefficients, mu, np.diag(scale**2), intercept, noise)


def gaussian_logistic(
    beta: ArrayLike, *, mean: ArrayLike | None = None, sigma: ArrayLike | float = 1.0, intercept: float = 0.0
) -> GaussianDGP:
    """Create Bernoulli(sigmoid(beta0+beta X)) DGP (Section C.1).

    Parameters
    ----------
    beta, mean, sigma : array_like
        Gaussian feature and logit parameters.
    intercept : float
        Logit intercept.

    Returns
    -------
    GaussianDGP
        Binary response model; analytic attribution targets logit output.
    """
    base = gaussian_linear(beta, mean=mean, sigma=sigma, intercept=intercept, noise=0)
    return GaussianDGP(base.beta, base.mean, base.covariance, intercept, 0, True)


def correlated_sem(beta: ArrayLike, rho: float = 0.5, *, noise: float = 1.0) -> GaussianDGP:
    """Create X2=rho X1+sqrt(1-rho^2) Z2 Gaussian SEM (Section C.1).

    assumption: remaining X are independent standard normals. Attribution
    remains marginal predictive attribution, not a causal-effect estimator.

    Parameters
    ----------
    beta : array_like of shape (p,), p >= 2
        Linear outcome coefficients.
    rho : float
        Correlation in [-1,1], including exact duplicate endpoints.
    noise : float
        Residual standard deviation.

    Returns
    -------
    GaussianDGP
        Analytically known correlated-feature DGP.
    """
    coefficients = finite_array(beta, "beta", 1)
    if len(coefficients) < 2 or not np.isfinite(rho) or abs(rho) > 1:
        raise ValueError("SEM requires p >= 2 and rho in [-1,1]")
    covariance = np.eye(len(coefficients))
    covariance[0, 1] = covariance[1, 0] = rho
    return GaussianDGP(coefficients, np.zeros_like(coefficients), covariance, noise=noise)


def exact_shapley(
    predict: Callable[[FloatArray], ArrayLike], x: ArrayLike, background: ArrayLike, *, max_features: int = 10
) -> FloatArray:
    """Enumerate marginal Shapley coalitions (Section 2.2 F / C.1).

    assumption: exact with respect to the supplied finite empirical background,
    not necessarily the population background. Suitable for small p only.

    Parameters
    ----------
    predict : callable
        Scalar-output model accepting a two-dimensional feature array.
    x, background : array_like
        Anchors (n,p) and finite reference sample (H,p).
    max_features : int
        Hard computational guard against accidental exponential execution.

    Returns
    -------
    ndarray of shape (n,p)
        Exact empirical-game signed attributions.
    """
    anchors = finite_array(x, "x", 2)
    bg = finite_array(background, "background", 2)
    n, p = anchors.shape
    if bg.shape[1] != p or p > max_features:
        raise ValueError("background feature mismatch or exact-enumeration limit exceeded")
    values = np.empty((1 << p, n))
    for mask in range(1 << p):
        columns = [j for j in range(p) if mask & (1 << j)]
        mixed = np.broadcast_to(bg, (n, *bg.shape)).copy()
        mixed[:, :, columns] = anchors[:, None, columns]
        outputs = finite_array(predict(mixed.reshape(-1, p)), "predict output", 1)
        if outputs.size != n * len(bg):
            raise ValueError("predict must return one scalar per input row")
        values[mask] = outputs.reshape(n, len(bg)).mean(axis=1)
    result = np.zeros((n, p))
    for j in range(p):
        for mask in range(1 << p):
            if not mask & (1 << j):
                size = mask.bit_count()
                weight = factorial(size) * factorial(p - size - 1) / factorial(p)
                result[:, j] += weight * (values[mask | (1 << j)] - values[mask])
    return result
