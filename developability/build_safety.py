"""CLI to compute ocular safety summary tables."""
from __future__ import annotations

import argparse
from pathlib import Path

from developability.faers import extract_ocular_events, load_openfda_events, summarise_ocular_events


OCULAR_SEARCH = "eye"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("ocular_safety.csv"))
    parser.add_argument("--term", default=OCULAR_SEARCH, help="Term to query from openFDA")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    events = load_openfda_events(args.term)
    ocular_events = extract_ocular_events(events)
    summary = summarise_ocular_events(ocular_events)
    summary.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
