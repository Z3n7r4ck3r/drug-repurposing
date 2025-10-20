"""Command-line utility to run propagation scoring for each disease."""
from __future__ import annotations

import argparse
import logging
import sqlite3
from pathlib import Path
from typing import Dict

import networkx as nx
import pandas as pd

from scoring.propagation import heat_diffusion


def load_graph(db_path: Path) -> nx.Graph:
    with sqlite3.connect(db_path) as conn:
        df = pd.read_sql_query("SELECT src, dst, sign FROM protein_edge", conn)
    graph = nx.DiGraph()
    for row in df.itertuples(index=False):
        weight = 1.0 if row.sign == "+" else -1.0
        graph.add_edge(row.src, row.dst, weight=weight)
    return graph


def load_seeds(path: Path) -> Dict[str, Dict[str, float]]:
    seeds_df = pd.read_csv(path)
    grouped = {}
    for disease_id, group in seeds_df.groupby("disease_id"):
        grouped[disease_id] = {row.gene_symbol: float(row.score) for row in group.itertuples()}
    return grouped


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, default=Path("data/knowledge_graph.sqlite"))
    parser.add_argument("--seeds", type=Path, default=Path("disease_seeds.csv"))
    parser.add_argument("--output", type=Path, default=Path("target_scores.csv"))
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    graph = load_graph(args.graph)
    logging.info("Loaded graph with %s nodes and %s edges", graph.number_of_nodes(), graph.number_of_edges())

    seed_sets = load_seeds(args.seeds)
    all_scores = []
    for disease_id, seeds in seed_sets.items():
        scores = heat_diffusion(graph, seeds)
        scores.insert(0, "disease_id", disease_id)
        all_scores.append(scores)

    if not all_scores:
        logging.warning("No scores generated")
        return

    result = pd.concat(all_scores, ignore_index=True)
    result.rename(columns={"node": "target_id"}, inplace=True)
    result.to_csv(args.output, index=False)
    logging.info("Wrote scores to %s", args.output)


if __name__ == "__main__":
    main()
