"""Network propagation utilities for disease module scoring."""
from __future__ import annotations

from typing import Mapping

import networkx as nx
import numpy as np
import pandas as pd


def heat_diffusion(
    graph: nx.Graph,
    seeds: Mapping[str, float],
    alpha: float = 0.7,
    tol: float = 1e-6,
    max_iter: int = 100,
) -> pd.DataFrame:
    """Perform heat diffusion propagation over the graph."""

    if graph.number_of_nodes() == 0:
        return pd.DataFrame(columns=["node", "score"])

    nodes = list(graph.nodes())
    index = {node: i for i, node in enumerate(nodes)}
    seed_vec = np.zeros(len(nodes))
    for node, weight in seeds.items():
        if node in index:
            seed_vec[index[node]] = weight

    adjacency = nx.to_numpy_array(graph, nodelist=nodes, weight="weight", dtype=float)
    degree = adjacency.sum(axis=1)
    degree[degree == 0] = 1
    normalised = adjacency / degree[:, None]

    scores = seed_vec.copy()
    for _ in range(max_iter):
        updated = alpha * seed_vec + (1 - alpha) * normalised.T.dot(scores)
        if np.linalg.norm(updated - scores, ord=1) < tol:
            scores = updated
            break
        scores = updated

    df = pd.DataFrame({"node": nodes, "score": scores})
    return df.sort_values("score", ascending=False).reset_index(drop=True)


def propagate_multiple(
    graph: nx.Graph,
    seed_sets: Mapping[str, Mapping[str, float]],
    **kwargs,
) -> pd.DataFrame:
    """Propagate multiple seed dictionaries and return stacked scores."""

    frames = []
    for label, seeds in seed_sets.items():
        result = heat_diffusion(graph, seeds, **kwargs)
        result.insert(0, "seed_set", label)
        frames.append(result)
    if not frames:
        return pd.DataFrame(columns=["seed_set", "node", "score"])
    return pd.concat(frames, ignore_index=True)
