"""Helpers for building differential expression signatures."""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence

import numpy as np
import pandas as pd


@dataclass
class ExpressionSignature:
    gene_symbol: str
    logfc: float
    pvalue: float

    def as_dict(self) -> dict:
        return {
            "gene_symbol": self.gene_symbol,
            "logfc": float(self.logfc),
            "pvalue": float(self.pvalue),
        }


def compute_logfc(
    expression: pd.DataFrame,
    case_samples: Sequence[str],
    control_samples: Sequence[str],
    gene_column: str = "gene_symbol",
) -> pd.DataFrame:
    """Compute log fold-changes between case and control cohorts."""

    missing = [col for col in list(case_samples) + list(control_samples) if col not in expression.columns]
    if missing:
        raise KeyError(f"Expression matrix missing columns: {missing}")

    if gene_column not in expression.columns:
        raise KeyError(f"Expression dataframe lacks '{gene_column}' column")

    matrix = expression.set_index(gene_column)
    case_mean = matrix[case_samples].astype(float).mean(axis=1)
    control_mean = matrix[control_samples].astype(float).mean(axis=1)
    logfc = np.log2((case_mean + 1e-6) / (control_mean + 1e-6))

    # Welch's t-test approximation using pandas/numpy to avoid scipy dependency
    case_var = matrix[case_samples].astype(float).var(axis=1, ddof=1)
    control_var = matrix[control_samples].astype(float).var(axis=1, ddof=1)
    case_n = len(case_samples)
    control_n = len(control_samples)

    se = np.sqrt(case_var / case_n + control_var / control_n)
    se.replace(0, np.nan, inplace=True)
    t_stat = logfc / se
    dof_num = (case_var / case_n + control_var / control_n) ** 2
    dof_den = ((case_var ** 2) / ((case_n ** 2) * (case_n - 1))) + (
        (control_var ** 2) / ((control_n ** 2) * (control_n - 1))
    )
    dof = dof_num / dof_den
    dof.replace(0, np.nan, inplace=True)

    # Survival function approximation for two-tailed p-value using Student's t distribution via scipy-like formula
    try:
        from scipy.stats import t as student_t  # type: ignore

        pvalues = 2 * student_t.sf(np.abs(t_stat), dof)
    except Exception:  # pragma: no cover - SciPy optional
        # Fallback using normal approximation
        from math import erf, sqrt

        pvalues = 2 * (1 - 0.5 * (1 + erf(np.abs(t_stat) / sqrt(2))))

    result = pd.DataFrame({
        gene_column: matrix.index,
        "logfc": logfc.values,
        "pvalue": pvalues,
    })
    return result.dropna()


def rank_signature(signature: pd.DataFrame, max_genes: int = 250) -> List[ExpressionSignature]:
    """Return a ranked list of expression signatures."""

    if signature.empty:
        return []

    ranked = signature.sort_values("logfc", ascending=False)
    top = ranked.head(max_genes)
    return [ExpressionSignature(row["gene_symbol"], row["logfc"], row["pvalue"]) for _, row in top.iterrows()]


def split_up_down(signature: pd.DataFrame, logfc_threshold: float = 1.0) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split the signature into up/down regulated subsets."""

    up = signature[signature["logfc"] >= logfc_threshold].copy()
    down = signature[signature["logfc"] <= -logfc_threshold].copy()
    return up, down
