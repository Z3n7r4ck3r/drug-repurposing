"""CLI helper to fetch and summarise ClinicalTrials.gov data."""
from __future__ import annotations

import argparse
from pathlib import Path

from developability.trials import fetch_trials, summarise_trials


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--condition", action="append", dest="conditions", required=True)
    parser.add_argument("--output", type=Path, default=Path("trial_summary.csv"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    trials = fetch_trials(args.conditions)
    summary = summarise_trials(trials)
    summary.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
