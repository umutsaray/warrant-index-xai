"""Section C.3 baselines and explicitly configured ERI toy stress protocol.

ERI source: https://arxiv.org/html/2602.05082v1, Definitions 5-9, Appendix G.
The aggregation rules are exact min/uniform mean. Synthetic stress transforms
below instantiate the paper's distances; they do NOT reproduce ERI-Bench or its
published numbers. WI components are never relabeled as ERI components.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from ._validation import FloatArray, finite_array, tolerance, unit_array
from .scenarios import Scenario

ERI_AXES = ("ERI-S", "ERI-R", "ERI-M", "ERI-D", "ERI-T")


def eri_min(axes: ArrayLike) -> FloatArray:
    """Aggregate supplied ERI axes by min (Section C.3; ERI Appendix G).

    Parameters
    ----------
    axes : array_like of shape (..., number_of_axes)
        ERI scores, not WI F/S/U/G/R. NaN is propagated.

    Returns
    -------
    ndarray
        Worst-axis score.
    """
    array = unit_array(axes, "ERI axes", allow_nan=True)
    if array.ndim < 1:
        raise ValueError("ERI axes need a final dimension")
    return array.min(axis=-1)


def eri_mean(axes: ArrayLike) -> FloatArray:
    """Aggregate supplied ERI axes by uniform mean (Section C.3; ERI Appendix G).

    Parameters
    ----------
    axes : array_like
        Predeclared ERI axes, including explicit missing axes if any.

    Returns
    -------
    ndarray
        Compensatory arithmetic mean without omitting missing scores.
    """
    array = unit_array(axes, "ERI axes", allow_nan=True)
    if array.ndim < 1:
        raise ValueError("ERI axes need a final dimension")
    return array.mean(axis=-1)


def single_component(z: ArrayLike, component: str) -> FloatArray:
    """Return F-only, S-only or G-only ablation (Section C.3).

    Parameters
    ----------
    z : array_like of shape (...,5)
        Actual WI component matrix.
    component : {'F','S','G'}
        Required ablation.

    Returns
    -------
    ndarray
        Selected unmodified score.
    """
    array = unit_array(z, "z", allow_nan=True)
    if array.ndim < 1 or array.shape[-1] != 5 or component not in {"F", "S", "G"}:
        raise ValueError("need five WI components and F, S or G ablation")
    return array[..., {"F": 0, "S": 1, "G": 3}[component]]


def eri_linear_stress(
    scenario: Scenario, *, scale: ArrayLike = 0.5, perturb_sd: float = 0.02, probes: int = 32, seed: int = 0
) -> dict[str, float]:
    """Instantiate five ERI axes for linear toy explanations (Section C.3).

    assumption: d(a,b)=RMS((a-b)/scale). S uses Gaussian input perturbations;
    R averages alpha in [.9,1) and collapses feature 1 into 0 in a linear
    model; M uses consecutive independent refits as permitted by ERI Def.7;
    D compares signed mean explanations under native populations with the
    SAME source model; T uses a synthetic smooth path, not real time series.
    Transform 1/(1+mean distance) precedes the outer averaging over examples.
    These probe budgets/distributions are configuration, not ERI-Bench defaults.

    Parameters
    ----------
    scenario : Scenario
        Linear coefficient/anchor fixture, including native population samples.
    scale : array_like
        Fixed positive distance scale per feature.
    perturb_sd : float
        Input perturbation standard deviation and smooth-path step scale.
    probes : int
        Number of Monte Carlo perturbation/redundancy/path steps.
    seed : int
        Independent probe seed.

    Returns
    -------
    dict
        Five scalar whole-explanation scores; not feature-wise ERI claims.
    """
    coefficients = finite_array(scenario.explanation_coefficients, "coefficients", 3)
    sites, t, p = coefficients.shape
    scales = tolerance(scale, p, "scale")
    if p < 2 or t < 2 or probes < 2 or not np.isfinite(perturb_sd) or perturb_sd < 0:
        raise ValueError("need p,T,probes >=2 and finite perturb_sd >=0")
    rng = np.random.default_rng(seed)
    delta = rng.normal(0, perturb_sd, size=(probes, p))
    s_dist = np.sqrt(np.mean((coefficients[:, :, None, :] * delta / scales) ** 2, axis=-1))
    s = float(np.mean(1 / (1 + s_dist.mean(axis=-1))))
    # Collapse E0+E1 vs collapsed linear model (c0+alpha*c1)*x0;
    # the discarded independent residual is the only difference.
    alpha = rng.uniform(0.9, 1, size=probes)
    residual = np.sqrt(1 - alpha**2) * rng.normal(size=probes)
    r_dist = np.abs(coefficients[:, :, None, 1] * residual) / scales[0] / np.sqrt(p - 1)
    r = float(np.mean(1 / (1 + r_dist.mean(axis=-1))))
    changes = np.diff(scenario.phi, axis=1)
    m_dist = np.sqrt(np.mean((changes / scales) ** 2, axis=-1))
    m = float(np.mean(1 / (1 + m_dist.mean(axis=1))))
    # ERI-D changes population with frozen source explainer; it does not
    # compare the independently fitted external model used by WI G.
    mean_x = scenario.native_anchors.mean(axis=1)
    d_dist = np.sqrt(np.mean((coefficients[0, None, :, :] * (mean_x[1:, None, :] - mean_x[0]) / scales) ** 2, axis=-1))
    d = float(np.mean(1 / (1 + d_dist)))
    path_step = rng.normal(0, perturb_sd, size=p)
    t_dist = np.sqrt(np.mean((coefficients * path_step / scales) ** 2, axis=-1))
    temporal = float(np.mean(1 / (1 + t_dist)))
    return dict(zip(ERI_AXES, [s, r, m, d, temporal], strict=True))
