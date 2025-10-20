"""Command-line utility to export harmonised drug-target mapping."""
from __future__ import annotations

import argparse
from pathlib import Path

from mapping.drug_target_map import attach_rxnorm_synonyms, export_mapping, load_graph_targets


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, default=Path("data/knowledge_graph.sqlite"))
    parser.add_argument("--output", type=Path, default=Path("drug_target_mapping.csv"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    targets = load_graph_targets(args.db)
    enriched = attach_rxnorm_synonyms(targets)
    export_mapping(enriched, args.output)


if __name__ == "__main__":
    main()
