"""Feature integration utilities for downstream modeling."""
from __future__ import annotations

import pandas as pd


def normalise_series(series: pd.Series) -> pd.Series:
    if series.empty:
        return series
    min_val = series.min()
    max_val = series.max()
    if min_val == max_val:
        return pd.Series([0.0] * len(series), index=series.index)
    return (series - min_val) / (max_val - min_val)


def integrate_features(
    target_scores: pd.DataFrame,
    drug_targets: pd.DataFrame,
    safety_signals: pd.DataFrame,
    trials: pd.DataFrame,
) -> pd.DataFrame:
    """Combine heterogeneous evidence tables into a single feature matrix."""

    features = target_scores.copy()
    if "node" in features.columns and "target_id" not in features.columns:
        features.rename(columns={"node": "target_id"}, inplace=True)
    if features.empty:
        features = pd.DataFrame(columns=["target_id", "propagation_score", "disease_id"])
    if "disease_id" not in features.columns:
        features["disease_id"] = ""
    if "score" in features.columns and "propagation_score" not in features.columns:
        features.rename(columns={"score": "propagation_score"}, inplace=True)

    if not drug_targets.empty:
        target_counts = (
            drug_targets.groupby("target_id")
            .agg(drug_count=("drug_id", "nunique"))
            .reset_index()
        )
        features = features.merge(target_counts, on="target_id", how="left")
    else:
        features["drug_count"] = 0

    if not safety_signals.empty and not drug_targets.empty:
        safety_scores = (
            safety_signals.groupby("drug")
            .agg(ocular_reports=("reports", "sum"))
            .reset_index()
        )
        target_map = (
            drug_targets[["drug_id", "target_id"]]
            .drop_duplicates()
            .merge(safety_scores, left_on="drug_id", right_on="drug", how="left")
        )
        target_safety = (
            target_map.groupby("target_id")
            .agg(ocular_reports=("ocular_reports", "sum"))
            .reset_index()
        )
        features = features.merge(target_safety, on="target_id", how="left")
        features["ocular_reports"].fillna(0, inplace=True)
    else:
        features["ocular_reports"] = 0

    if not trials.empty:
        trial_counts = (
            trials.assign(condition=trials["condition"] if "condition" in trials.columns else trials.get("conditions", ""))
        )
        trial_counts["active_trials"] = trial_counts.get("active_trials", trial_counts.get("nct_id", 0))
        trial_counts = trial_counts[["condition", "active_trials"]].drop_duplicates()
        features = features.merge(
            trial_counts.rename(columns={"condition": "disease_id"}),
            on="disease_id",
            how="left",
        )
        features["active_trials"].fillna(0, inplace=True)
    else:
        features["active_trials"] = 0

    for col in ["propagation_score", "drug_count", "ocular_reports", "active_trials"]:
        if col in features.columns:
            features[f"{col}_norm"] = normalise_series(features[col].fillna(0))

    return features
