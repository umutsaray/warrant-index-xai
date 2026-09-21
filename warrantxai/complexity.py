"""Section E: finite-population subsampling and content-addressed caching."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re

import numpy as np
from numpy.typing import ArrayLike

from ._validation import FloatArray, finite_array, tolerance


def sample_pairs(size: int, count: int, seed: int | np.random.Generator = 0) -> np.ndarray:
    """Sample unordered pairs uniformly without replacement (Section E.2).

    Parameters
    ----------
    size : int
        Number of objects, at least two.
    count : int
        Between one and size*(size-1)/2 inclusive.
    seed : int or Generator
        Randomness source.

    Returns
    -------
    ndarray of shape (count,2)
        Unique pairs i<j, without constructing the quadratic pair list.
    """
    total = size * (size - 1) // 2
    if size < 2 or count < 1 or count > total:
        raise ValueError("invalid pair population/count")
    rng = np.random.default_rng(seed)
    ranks = rng.choice(total, count, replace=False)
    ends = np.cumsum(np.arange(size - 1, 0, -1, dtype=np.int64))
    first = np.searchsorted(ends, ranks, side="right")
    starts = np.concatenate(([0], ends[:-1]))
    second = first + 1 + ranks - starts[first]
    return np.column_stack([first, second])


def incomplete_stability(
    phi: ArrayLike, a: ArrayLike, training_pairs: int, seed: int | np.random.Generator = 0
) -> FloatArray:
    """Approximate S with uniformly sampled training pairs (Section E.2).

    assumption: conditional on the cached tensor, site-level means are unbiased
    finite-population estimates. The minimum across sites is NOT unbiased.
    This estimates the unchanged Section 2.2 formula, not a replacement score.

    Parameters
    ----------
    phi : array_like of shape (sites,T,n,p)
        Cached mean attributions.
    a : array_like
        Fixed feature scales.
    training_pairs : int
        Pair budget; all pairs recovers exact S.
    seed : int or Generator
        Independent sampling source.

    Returns
    -------
    ndarray of shape (p,)
        Approximate worst-site S.
    """
    array = finite_array(phi, "phi", 4)
    scales = tolerance(a, array.shape[-1], "a")
    pairs = sample_pairs(array.shape[1], training_pairs, seed)
    loss = np.zeros((array.shape[0], array.shape[-1]))
    for first, second in pairs:
        loss += (np.minimum(1, np.abs((array[:, first] - array[:, second]) / scales)) ** 2).mean(axis=1)
    return (1 - loss / len(pairs)).min(axis=0)


def incomplete_rank(
    imp: ArrayLike, eta: float, training_pairs: int, comparisons_per_feature: int, seed: int | np.random.Generator = 0
) -> FloatArray:
    """Approximate feature-wise R with paired subsampling (Section E.2).

    Every feature samples its own uniformly chosen other features. Each is
    represented; top-k filtering would change the target and is not used.
    Sitewise means are unbiased conditional on the tensor, their min is not.

    Parameters
    ----------
    imp : array_like of shape (sites,T,p)
        Cached mean absolute importance.
    eta : float
        Fixed nonnegative tie tolerance.
    training_pairs, comparisons_per_feature : int
        Sampling budgets; full budgets recover exact R.
    seed : int or Generator
        Randomness source.

    Returns
    -------
    ndarray of shape (p,)
        Approximate R for all original features.
    """
    array = finite_array(imp, "importance", 3)
    p = array.shape[-1]
    if np.any(array < 0) or not np.isfinite(eta) or eta < 0 or not 1 <= comparisons_per_feature <= p - 1:
        raise ValueError("need nonnegative importance/eta and 1<=comparisons<=p-1")
    rng = np.random.default_rng(seed)
    pairs = sample_pairs(array.shape[1], training_pairs, rng)
    neighbors = np.stack([rng.choice(p - 1, comparisons_per_feature, replace=False) for _ in range(p)])
    neighbors += neighbors >= np.arange(p)[:, None]
    difference = array[:, :, :, None] - array[:, :, neighbors]
    order = (difference > eta).astype(np.int8) - (difference < -eta).astype(np.int8)
    loss = np.zeros((array.shape[0], p))
    for first, second in pairs:
        loss += (np.abs(order[:, first] - order[:, second]) / 2).mean(axis=-1)
    return (1 - loss / len(pairs)).min(axis=0)


def stratified_anchor_subsample(
    labels: ArrayLike, count: int, seed: int | np.random.Generator = 0
) -> tuple[np.ndarray, FloatArray]:
    """Subsample anchors while retaining every stratum (Section E.2).

    assumption: stratum quotas are chosen by population counts, never scores.
    Returned weights N_h/(N*m_h) estimate population means; simply averaging
    the selected rows can change the estimand when allocation is disproportionate.

    Parameters
    ----------
    labels : array_like of shape (n,)
        Complete finite stratum labels (strings or numbers).
    count : int
        Sample count >= number of strata and <= n.
    seed : int or Generator
        Sampling source.

    Returns
    -------
    indices, mean_weights : ndarray
        Selected row indices and population-mean weights summing to one.
    """
    strata = np.asarray(labels)
    if strata.ndim != 1 or not strata.size:
        raise ValueError("labels must be nonempty vector")
    if strata.dtype.kind in "fc" and not np.isfinite(strata).all():
        raise ValueError("labels must not be missing/infinite")
    _, inverse, population = np.unique(strata, return_inverse=True, return_counts=True)
    if not len(population) <= count <= len(strata):
        raise ValueError("count must cover each stratum without exceeding n")
    ideal = count * population / population.sum()
    quotas = np.ones(len(population), dtype=int)
    for _ in range(count - len(population)):
        score = np.where(quotas < population, ideal - quotas, -np.inf)
        quotas[np.argmax(score)] += 1
    rng = np.random.default_rng(seed)
    indices = np.concatenate(
        [rng.choice(np.flatnonzero(inverse == h), quotas[h], replace=False) for h in range(len(population))]
    )
    weights = np.concatenate(
        [np.full(quotas[h], population[h] / (len(strata) * quotas[h])) for h in range(len(population))]
    )
    return indices, weights


class AttributionCache:
    """Content-addressed numeric attribution cache (Section E.2).

    Keys include anchors/background/model arrays and explicit semantic metadata.
    No pickle loading. Atomic replacement prevents partial cached tensors.
    """

    def __init__(self, directory: str | Path) -> None:
        """Configure a local cache directory (Section E.2).

        Parameters
        ----------
        directory : str or Path
            Dedicated cache location, created on first save.
        """
        self.directory = Path(directory)

    @staticmethod
    def key(*arrays: ArrayLike, metadata: dict) -> str:
        """Hash model/anchor/background content plus semantics (Section E.2).

        Parameters
        ----------
        arrays : array_like
            All numerical inputs affecting attributions, including coefficients.
        metadata : dict
            JSON-serializable protocol/explainer/output-scale/version settings.

        Returns
        -------
        str
            SHA256 cache identity; callers must include every semantic dependency.
        """
        digest = hashlib.sha256(json.dumps(metadata, sort_keys=True, allow_nan=False).encode())
        for value in arrays:
            array = np.ascontiguousarray(finite_array(value, "cache input"))
            digest.update(str(array.shape).encode())
            digest.update(array.dtype.str.encode())
            digest.update(array.tobytes())
        return digest.hexdigest()

    def _path(self, key: str) -> Path:
        """Reject invalid/traversing keys (Section E.2).

        Returns
        -------
        Path
            Valid cache-file location.
        """
        if not re.fullmatch(r"[0-9a-f]{64}", key):
            raise ValueError("cache key must be SHA256 hexadecimal")
        return self.directory / f"{key}.npy"

    def get(self, key: str) -> FloatArray | None:
        """Read a cached attribution tensor without pickle (Section E.2).

        Returns
        -------
        ndarray or None
            Validated numeric data, or cache miss.
        """
        path = self._path(key)
        if not path.exists():
            return None
        return finite_array(np.load(path, allow_pickle=False), "cached attributions")

    def put(self, key: str, values: ArrayLike) -> None:
        """Atomically write numeric cache content (Section E.2).

        Parameters
        ----------
        key : str
            SHA256 identity.
        values : array_like
            Finite attribution tensor.

        Returns
        -------
        None
        """
        import tempfile

        path = self._path(key)
        array = finite_array(values, "attributions")
        self.directory.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=self.directory, suffix=".npy", delete=False) as handle:
            temporary = Path(handle.name)
            np.save(handle, array, allow_pickle=False)
        try:
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)
