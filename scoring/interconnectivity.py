"""Graph interconnectivity metrics used as complementary evidence."""
from __future__ import annotations

from itertools import combinations
from typing import Iterable

import networkx as nx
import numpy as np
import pandas as pd


def average_shortest_path(graph: nx.Graph, nodes: Iterable[str]) -> float:
    sub_nodes = [node for node in nodes if node in graph]
    if len(sub_nodes) < 2:
        return float("nan")
    lengths = []
    for a, b in combinations(sub_nodes, 2):
        try:
            lengths.append(nx.shortest_path_length(graph, a, b, weight="weight"))
        except nx.NetworkXNoPath:
            continue
    if not lengths:
        return float("inf")
    return float(np.mean(lengths))


def connectivity_zscore(graph: nx.Graph, nodes: Iterable[str], iterations: int = 1000) -> float:
    """Compare the observed average path length against random node sets."""

    observed = average_shortest_path(graph, nodes)
    if np.isnan(observed) or np.isinf(observed):
        return float("nan")

    sub_nodes = [node for node in nodes if node in graph]
    if len(sub_nodes) < 2:
        return float("nan")

    nodes_list = list(graph.nodes())
    rng = np.random.default_rng(42)
    random_lengths = []
    for _ in range(iterations):
        sample = rng.choice(nodes_list, size=len(sub_nodes), replace=False)
        random_lengths.append(average_shortest_path(graph, sample))

    random_array = np.array(random_lengths, dtype=float)
    mu = np.nanmean(random_array)
    sigma = np.nanstd(random_array)
    if sigma == 0 or np.isnan(sigma):
        return float("nan")
    return float((observed - mu) / sigma)


def enrich_interactions(graph: nx.Graph, nodes: Iterable[str]) -> pd.DataFrame:
    """Return edges among the provided nodes for inspection."""

    sub_nodes = [node for node in nodes if node in graph]
    subgraph = graph.subgraph(sub_nodes)
    records = []
    for u, v, data in subgraph.edges(data=True):
        records.append({"src": u, "dst": v, "weight": data.get("weight", 1.0)})
    return pd.DataFrame(records)
