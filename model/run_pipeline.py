"""End-to-end model feature integration and optional training."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from model.integrate_features import integrate_features
from model.train_validate import TrainingResult, train_model


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=Path("target_scores.csv"))
    parser.add_argument("--mapping", type=Path, default=Path("drug_target_mapping.csv"))
    parser.add_argument("--safety", type=Path, default=Path("ocular_safety.csv"))
    parser.add_argument("--trials", type=Path, default=Path("trial_summary.csv"))
    parser.add_argument("--labels", type=Path, help="Optional training labels CSV")
    parser.add_argument("--features-out", type=Path, default=Path("model_features.csv"))
    parser.add_argument("--metrics-out", type=Path, default=Path("model_metrics.json"))
    return parser.parse_args()


def load_csv(path: Path) -> pd.DataFrame:
    if path is None or not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def main() -> None:
    args = parse_args()

    scores = load_csv(args.scores)
    mapping = load_csv(args.mapping)
    safety = load_csv(args.safety)
    trials = load_csv(args.trials)

    features = integrate_features(scores, mapping, safety, trials)

    if args.labels and args.labels.exists():
        labels = pd.read_csv(args.labels)
        features = features.merge(labels, on=["target_id", "disease_id"], how="left")
        features["label"].fillna(0, inplace=True)

        result: TrainingResult = train_model(features)
        args.metrics_out.write_text(json.dumps({"auc": result.auc, "average_precision": result.average_precision}))

    features.to_csv(args.features_out, index=False)


if __name__ == "__main__":
    main()
