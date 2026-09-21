"""Warrant Index: fixed-formula, synthetic-first explanation audits."""

from .aggregation import adequate, warrant_index, claim_index
from .components import component_scores
from .api import WarrantIndexEvaluator

__version__ = "0.1.0"
__all__ = ["adequate", "warrant_index", "claim_index", "component_scores", "WarrantIndexEvaluator"]
