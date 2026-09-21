"""Section 2/Phase 2 API: sklearn-style estimator for fitted linear models.

The requested fit(models, data_by_site) contract is an audit API, not the usual
supervised fit(X,y) contract. It supports sklearn get_params/set_params/clone,
but does NOT claim drop-in supervised Pipeline/GridSearchCV compatibility.
"""

from __future__ import annotations

from typing import Any, Callable, Sequence

import numpy as np
from numpy.typing import ArrayLike
import pandas as pd
from sklearn.base import BaseEstimator
from sklearn.utils.validation import check_is_fitted

from ._validation import finite_array, tolerance
from .aggregation import adequate, claim_index, warrant_index
from .components import COMPONENT_NAMES, component_scores, importance, importance_ci, pairwise_orders


class WarrantIndexEvaluator(BaseEstimator):
    """Audit fitted linear/logit model ensembles (Section 2.2-2.3).

    Parameters
    ----------
    a, b, u : array_like
        Fixed positive tolerances. Never inferred from evaluation data.
    eta : float
        Fixed rank tie tolerance.
    tau : array_like
        Scalar or five strict adequacy thresholds.
    background_mean : array_like, optional
        Calibrated background mean; assumption: zero vector when omitted.
    anchors : array_like, optional
        Common Q anchors; assumption: source-site data used if omitted.
    common_support : bool or array_like of bool
        Caller assertion about every external site; no automatic overlap test.
    alpha_U : float
        Conditional normal-CI error level.
    claim_features : sequence of int, optional
        Preregistered J*. None reports feature scores without a claim score.
    attributor : callable, optional
        (model,anchors,background_mean)->(n,p) mean attribution. Default closed
        form. Reference always uses the model's exact marginal linear/logit game.
    generalizability_mode : {"claim", "magnitude_only"}
        Default signed-and-magnitude claim; magnitude_only is a sensitivity.
    """

    def __init__(
        self,
        *,
        a: ArrayLike = 0.5,
        b: ArrayLike = 0.5,
        u: ArrayLike = 3.0,
        eta: float = 0.02,
        tau: ArrayLike = 0.8,
        background_mean: ArrayLike | None = None,
        anchors: ArrayLike | None = None,
        common_support: ArrayLike = True,
        alpha_U: float = 0.05,
        claim_features: Sequence[int] | None = None,
        attributor: Callable | None = None,
        generalizability_mode: str = "claim",
    ) -> None:
        """Store constructor parameters unchanged for sklearn cloning (Section 2).

        Returns
        -------
        None
        """
        self.a = a
        self.b = b
        self.u = u
        self.eta = eta
        self.tau = tau
        self.background_mean = background_mean
        self.anchors = anchors
        self.common_support = common_support
        self.alpha_U = alpha_U
        self.claim_features = claim_features
        self.attributor = attributor
        self.generalizability_mode = generalizability_mode

    def fit(self, models: Sequence[Sequence[Any]], data_by_site: Sequence[ArrayLike]) -> WarrantIndexEvaluator:
        """Cache explanations of ALREADY fitted model replicas (Section 2.2).

        assumption: model coefficients describe a scalar linear predictor, or a
        binary logit. Nonlinear feature preprocessing and probability-output
        logistic explanations are outside this adapter's contract. Supplied
        site data must be independent audit data, not reused training samples.
        Site-specific datasets are validated; scoring uses one shared anchor Q.

        Parameters
        ----------
        models : sequence of sequences
            Equal-length fitted model lists, one per site; T >= 2. Each has
            coef_ of shape (p,) or (1,p). No model is trained or modified here.
        data_by_site : sequence of array_like
            Each site's evaluation data (n_d,p). Site zero is source.

        Returns
        -------
        WarrantIndexEvaluator
            Fitted evaluator with attributions and exact linear references.
        """
        # Invalidate the previous audit before attempting a new fit.
        for name in list(vars(self)):
            if name.endswith("_"):
                delattr(self, name)
        if len(models) == 0 or len(models) != len(data_by_site):
            raise ValueError("models and data_by_site must have matching nonzero site counts")
        data = [finite_array(x, "site data", 2).copy() for x in data_by_site]
        p = data[0].shape[1]
        if p < 2 or any(x.shape[1] != p for x in data):
            raise ValueError("aligned p >= 2 features required at all sites")
        t = len(models[0])
        if t < 2 or any(len(group) != t for group in models):
            raise ValueError("equal T >= 2 fitted replicas required at every site")
        for name in ["a", "b", "u"]:
            tolerance(getattr(self, name), p, name)
        background = (
            np.zeros(p) if self.background_mean is None else finite_array(self.background_mean, "background_mean", 1)
        )
        anchors = data[0] if self.anchors is None else finite_array(self.anchors, "anchors", 2)
        if background.shape != (p,) or anchors.shape[1] != p:
            raise ValueError("background/anchor feature mismatch")
        reference = np.empty((len(models), t, len(anchors), p))
        candidate = np.empty_like(reference)
        for d, group in enumerate(models):
            for rep, model in enumerate(group):
                if not hasattr(model, "coef_"):
                    raise ValueError("adapter needs fitted scalar linear/logit models with coef_")
                coefficients = finite_array(model.coef_, "model.coef_")
                if coefficients.shape == (1, p):
                    coefficients = coefficients[0]
                if coefficients.shape != (p,):
                    raise ValueError("multiclass/multioutput coefficients are not supported")
                reference[d, rep] = (anchors - background) * coefficients
                value = (
                    reference[d, rep]
                    if self.attributor is None
                    else self.attributor(model, anchors.copy(), background.copy())
                )
                attribution = finite_array(value, "attributor output", 2)
                if attribution.shape != (len(anchors), p):
                    raise ValueError("attributor must return (n,p)")
                candidate[d, rep] = attribution
        lo, hi = importance_ci(candidate, self.alpha_U)
        z = component_scores(
            candidate,
            reference,
            lo,
            hi,
            a=self.a,
            b=self.b,
            u=self.u,
            eta=self.eta,
            common_support=self.common_support,
            generalizability_mode=self.generalizability_mode,
        )
        # Validate threshold and claim contract before recording successful fit.
        adequate(z, self.tau)
        if self.claim_features is not None:
            claim_index(warrant_index(z), self.claim_features)
        self.n_features_in_ = p
        self.attributions_ = candidate
        self.reference_ = reference
        self.lower_ = lo
        self.upper_ = hi
        self.anchors_ = anchors.copy()
        self.background_mean_ = background.copy()
        self.site_sizes_ = tuple(len(x) for x in data)
        self.components_ = z
        self.generalizability_mode_ = self.generalizability_mode
        self.component_names_ = (
            *COMPONENT_NAMES[:3],
            "G_claim" if self.generalizability_mode == "claim" else "G_magnitude_only_sensitivity",
            COMPONENT_NAMES[4],
        )
        return self

    def evaluate(self) -> dict[str, Any]:
        """Return point scores and explicit assessment status (Section 2.3).

        Returns
        -------
        dict
            Components, WI, nullable point adequacy, status and optional
            WI_claim. No automatic inferential certification is performed.
        """
        check_is_fitted(self, "components_")
        wi = warrant_index(self.components_)
        return {
            "components": self.components_.copy(),
            "component_names": self.component_names_,
            "generalizability_mode": self.generalizability_mode_,
            "WI": wi,
            "point_adequate": adequate(self.components_, self.tau),
            "status": np.where(np.isfinite(wi), "assessed", "not assessable"),
            "WI_claim": None if self.claim_features is None else claim_index(wi, self.claim_features),
        }

    def feature_report(self) -> pd.DataFrame:
        """Return feature scores plus non-informative-tie diagnostics (Section 2).

        Returns
        -------
        DataFrame
            F/S/U/G/R, WI, source importance, point adequacy, missingness,
            bottleneck, CI summaries and non-tied pair fraction.
        """
        result = self.evaluate()
        z = result["components"]
        imp = importance(self.attributions_)
        orders = pairwise_orders(imp, self.eta)
        frame = pd.DataFrame(z, columns=self.component_names_)
        frame.index.name = "feature"
        frame["source_importance"] = imp[0].mean(0)
        frame["WI"] = result["WI"]
        frame["point_adequate"] = pd.array(
            [pd.NA if np.isnan(x) else bool(x) for x in result["point_adequate"]], dtype="boolean"
        )
        frame["status"] = result["status"]
        frame["bottleneck"] = [
            "not assessable" if np.isnan(row).any() else self.component_names_[int(np.argmin(row))] for row in z
        ]
        frame["source_CI_lower_mean"] = self.lower_[0].mean(0)
        frame["source_CI_upper_mean"] = self.upper_[0].mean(0)
        frame["non_tied_fraction"] = np.mean((orders != 0).sum(-1) / (self.n_features_in_ - 1), axis=(0, 1))
        frame["certification"] = "not performed: supply valid component p-values"
        return frame
