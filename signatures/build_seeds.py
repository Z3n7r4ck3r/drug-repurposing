"""Command-line helper to build disease seed tables."""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from signatures.seed_builder import (
    SeedEvidence,
    assemble_seed_table,
    load_gwas_seeds,
    load_opentargets_seeds,
    summarise_by_disease,
    save_seed_table,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("disease_seeds.csv"))
    parser.add_argument("--summary", type=Path, default=Path("disease_seed_summary.csv"))
    parser.add_argument(
        "--gwas", action="append", dest="gwas_ids", default=[], help="EFO identifiers to query from GWAS Catalog"
    )
    parser.add_argument("--min-score", type=float, default=0.2, help="Minimum OpenTargets score")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO)

    opentargets = load_opentargets_seeds(min_score=args.min_score)
    gwas = load_gwas_seeds(args.gwas_ids) if args.gwas_ids else None

    curated = [
        SeedEvidence("MONDO:0005015", "CFH", 1.0, "Curated", "literature", "AMD seed")
    ]

    seeds = assemble_seed_table(opentargets, gwas, curated)
    save_seed_table(seeds, str(args.output))

    summary = summarise_by_disease(seeds)
    summary.to_csv(args.summary, index=False)
    logging.info("Wrote summary to %s", args.summary)


if __name__ == "__main__":
    main()
