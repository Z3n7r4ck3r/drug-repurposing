import math

import pytest

pd = pytest.importorskip("pandas")

from signatures.disease_module import (
    ExpressionStudy,
    integrate_disease_evidence,
    meta_analyse_studies,
    run_study_analysis,
)


def _mock_expression(study_id: str, case_values, control_values) -> ExpressionStudy:
    genes = ["GENE1", "GENE2", "GENE3"]
    case_cols = [f"{study_id}_case_{idx}" for idx in range(len(case_values[0]))]
    control_cols = [f"{study_id}_control_{idx}" for idx in range(len(control_values[0]))]
    data = {"gene_symbol": genes}
    for idx, col in enumerate(case_cols):
        data[col] = [vals[idx] for vals in case_values]
    for idx, col in enumerate(control_cols):
        data[col] = [vals[idx] for vals in control_values]
    expression = pd.DataFrame(data)
    return ExpressionStudy(
        study_id=study_id,
        disease_id="D1",
        expression=expression,
        case_samples=case_cols,
        control_samples=control_cols,
    )


def test_run_study_analysis_produces_effect_sizes():
    study = _mock_expression("S1", case_values=[[5, 6], [3, 3], [2, 2]], control_values=[[2, 2], [2, 2], [2, 2]])
    tables = run_study_analysis([study])
    assert len(tables) == 1
    result = tables[0]
    assert set(result.columns) == {"disease_id", "study_id", "gene_symbol", "logfc", "se", "pvalue"}
    assert math.isclose(result.loc[result["gene_symbol"] == "GENE1", "logfc"].iloc[0], math.log2(5.5 / 2.0), rel_tol=1e-4)


def test_meta_analyse_combines_two_studies():
    study1 = _mock_expression("S1", case_values=[[5, 6], [3, 3], [2, 2]], control_values=[[2, 2], [2, 2], [2, 2]])
    study2 = _mock_expression("S2", case_values=[[6, 6], [2, 2], [2, 2]], control_values=[[2, 2], [2, 2], [2, 2]])
    tables = run_study_analysis([study1, study2])
    meta = meta_analyse_studies(tables)
    gene1 = meta.loc[(meta["disease_id"] == "D1") & (meta["gene_symbol"] == "GENE1")]
    assert not gene1.empty
    assert gene1["n_studies"].iloc[0] == 2
    assert math.isfinite(gene1["meta_logfc"].iloc[0])
    assert 0 <= gene1["meta_pvalue"].iloc[0] <= 1


def test_integrate_disease_evidence_merges_expression_and_seeds():
    meta = pd.DataFrame(
        {
            "disease_id": ["D1", "D1"],
            "gene_symbol": ["GENE1", "GENE2"],
            "meta_logfc": [1.2, -0.5],
            "meta_pvalue": [1e-4, 1e-2],
            "n_studies": [2, 2],
        }
    )
    seeds = pd.DataFrame(
        {
            "disease_id": ["D1"],
            "gene_symbol": ["GENE2"],
            "score": [0.8],
        }
    )
    integrated = integrate_disease_evidence(meta, seeds)
    assert set(integrated.columns) == {
        "disease_id",
        "gene_symbol",
        "combined_score",
        "expression_score",
        "seed_score",
        "n_studies",
        "sources",
    }
    gene2 = integrated[integrated["gene_symbol"] == "GENE2"].iloc[0]
    assert gene2["seed_score"] == 0.8
    assert "seed" in gene2["sources"]
    assert 0 <= gene2["combined_score"] <= 1
