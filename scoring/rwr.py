"""Random walk with restart implementation for disease module scoring."""
from __future__ import annotations

from typing import Mapping

import networkx as nx
import numpy as np
import pandas as pd


def random_walk_with_restart(
    graph: nx.Graph,
    seeds: Mapping[str, float],
    restart_prob: float = 0.3,
    tol: float = 1e-6,
    max_iter: int = 200,
) -> pd.DataFrame:
    if graph.number_of_nodes() == 0:
        return pd.DataFrame(columns=["node", "score"])

    nodes = list(graph.nodes())
    index = {node: i for i, node in enumerate(nodes)}
    seed_vec = np.zeros(len(nodes))
    norm = sum(seeds.values()) or 1.0
    for node, weight in seeds.items():
        if node in index:
            seed_vec[index[node]] = weight / norm

    adjacency = nx.to_numpy_array(graph, nodelist=nodes, weight="weight", dtype=float)
    degree = adjacency.sum(axis=1)
    degree[degree == 0] = 1
    transition = adjacency / degree[:, None]

    scores = seed_vec.copy()
    for _ in range(max_iter):
        updated = (1 - restart_prob) * transition.T.dot(scores) + restart_prob * seed_vec
        if np.linalg.norm(updated - scores, ord=1) < tol:
            scores = updated
            break
        scores = updated

    df = pd.DataFrame({"node": nodes, "score": scores})
    return df.sort_values("score", ascending=False).reset_index(drop=True)
