
"""rnanorm.deseq2

Minimal DESeq2 *median‑of‑ratios* normalisation transformer for RNAnorm.

The implementation mirrors the Bioconductor reference in a handful of clear
steps:

1. **Replace zero counts with *NaN*** so they are excluded from *all* further
   calculations (this reproduces DESeq2's default *ratio* method that drops
   genes containing even a single zero).
2. **Learn the pseudo‑sample** – per‑gene *log* geometric means calculated with
   ``np.mean`` (not ``np.nanmean``!) so that the presence of *any* *NaN*
   propagates and the whole gene is discarded.
3. **Sample‑specific size factor** – ``sf = exp(nanmedian(log(counts) − logG))``.
4. **Normalised counts** – divide raw counts by the size factor (no extra
   rescaling because downstream ML models expect the original DESeq2 scale).

Example
~~~~~~~
>>> from rnanorm.datasets import load_toy_data
>>> from rnanorm import DESeq2
>>> X = load_toy_data().exp
>>> de = DESeq2().set_output(transform="pandas")
>>> X_train, X_test = X.iloc[:3], X.iloc[3:]
>>> X_train_norm = de.fit_transform(X_train)
>>> X_test_norm  = de.transform(X_test)
"""

# from __future__ import annotations was previously scheduled to become mandatory in Python 3.10,
from __future__ import annotations

import typing
from typing import Any, Callable, Optional
import numpy as np
import pandas as pd
import numpy.typing as npt
from pandera.typing import DataFrame, Series
from sklearn.base import BaseEstimator, OneToOneFeatureMixin, TransformerMixin
from sklearn.utils._set_output import _get_output_config
from sklearn.utils.validation import check_is_fitted, validate_data

# the code in typing --------------------
#: Type for a 1D matrix of numeric values.
Numeric1D = typing.Union[Series[float], npt.NDArray[np.float64]]
#: Type for a 2D matrix of numeric values.
Numeric2D = typing.Union[DataFrame[float], npt.NDArray[np.float64]]

if hasattr(typing, "Self"):
    Self = typing.Self
else:
    Self = object

class DESeq2(OneToOneFeatureMixin, TransformerMixin, BaseEstimator):
    """DESeq2 *median‑of‑ratios* normalisation.

    Parameters
    ----------
    locfunc
        Location estimator applied over genes (axis = 1).  Must accept an
        ``axis`` argument like the NumPy reduction APIs.  The default
        ``numpy.nanmedian`` faithfully reproduces DESeq2.
    dtype
        Floating‑point dtype used internally (defaults to ``float64``).
    """

    def __init__(self, *, locfunc: Callable[..., np.ndarray] = np.nanmedian, dtype: np.dtype = np.float64) -> None:
        self.locfunc = locfunc
        self.dtype = dtype

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    @staticmethod
    def _log_counts(X: np.ndarray, *, dtype: np.dtype) -> np.ndarray:
        """Return natural‑log counts with zeros turned into *NaN* (vectorised)."""
        with np.errstate(divide="ignore"):
            return np.where(X > 0, np.log(X, dtype=dtype), np.nan)

    # ------------------------------------------------------------------
    # scikit‑learn API
    # ------------------------------------------------------------------

    def fit(self, X: Numeric2D, y: Optional[Numeric1D] = None, **fit_params: Any) -> "DESeq2":
        """Learn gene‑wise log‑geometric means from the training set.

        Genes that contain *any* zero become entirely *NaN* and are excluded
        from downstream size‑factor calculations – just like the reference
        DESeq2 implementation.
        """
        X = validate_data(self, X, ensure_all_finite=False, reset=True, dtype=self.dtype)
        logX = self._log_counts(X, dtype=self.dtype)
        # Use np.mean (not nanmean) so that a single NaN propagates.
        self.loggeomeans_ = np.mean(logX, axis=0, dtype=self.dtype)
        return self

    def get_size_factors(self, X: Numeric2D) -> Numeric1D:
        """Compute size factors for *X* without modifying the counts."""
        check_is_fitted(self, ["loggeomeans_"])
        X = validate_data(self, X, ensure_all_finite=False, reset=False, dtype=self.dtype)
        logX = self._log_counts(X, dtype=self.dtype)
        diff = logX - self.loggeomeans_  # broadcast over samples
        sf = np.exp(self.locfunc(diff, axis=1)).astype(self.dtype)

        cfg = _get_output_config("transform", self)
        if cfg.get("dense") == "pandas" and isinstance(X, pd.DataFrame):
            return pd.Series(sf, index=X.index)
        return sf

    def transform(self, X: Numeric2D) -> Numeric2D:
        """Return DESeq2‑normalised counts for *X*."""
        sf = self.get_size_factors(X)
        if isinstance(sf, pd.Series):
            sf = sf.to_numpy()
        X = validate_data(self, X, ensure_all_finite=False, reset=False, dtype=self.dtype)
        normed = X / sf[:, None]

        cfg = _get_output_config("transform", self)
        if cfg.get("dense") == "pandas" and isinstance(X, pd.DataFrame):
            return pd.DataFrame(normed, index=X.index, columns=X.columns)
        return normed
