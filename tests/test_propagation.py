import networkx as nx

from scoring.propagation import heat_diffusion


def test_heat_diffusion_returns_scores():
    graph = nx.Graph()
    graph.add_edge("A", "B", weight=1.0)
    graph.add_edge("B", "C", weight=1.0)

    scores = heat_diffusion(graph, {"A": 1.0})
    assert not scores.empty
    assert scores.loc[scores["node"] == "A", "score"].iloc[0] >= scores["score"].min()
