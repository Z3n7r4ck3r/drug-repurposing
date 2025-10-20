"""Training utilities for the repurposing ranking model."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import train_test_split


@dataclass
class TrainingResult:
    model: LogisticRegression
    auc: float
    average_precision: float


FEATURE_COLUMNS = [
    "propagation_score_norm",
    "drug_count_norm",
    "ocular_reports_norm",
    "active_trials_norm",
]


def prepare_dataset(features: pd.DataFrame, label_column: str = "label") -> Tuple[pd.DataFrame, pd.Series]:
    if label_column not in features.columns:
        raise KeyError(f"Label column '{label_column}' missing from features")
    X = features[[col for col in FEATURE_COLUMNS if col in features.columns]].fillna(0)
    y = features[label_column].astype(int)
    return X, y


def train_model(
    features: pd.DataFrame,
    label_column: str = "label",
    test_size: float = 0.2,
    random_state: int = 42,
) -> TrainingResult:
    X, y = prepare_dataset(features, label_column)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_state)

    model = LogisticRegression(max_iter=1000)
    model.fit(X_train, y_train)

    probs = model.predict_proba(X_test)[:, 1]
    auc = roc_auc_score(y_test, probs)
    ap = average_precision_score(y_test, probs)
    return TrainingResult(model=model, auc=auc, average_precision=ap)
