from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

pytest.importorskip("pandas")
import pandas as pd

from graph import build_graph


def _write_dataframe(path: Path, df: pd.DataFrame) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path


def _write_tsv(path: Path, df: pd.DataFrame) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    return path


def test_assemble_tables_with_local_sources(monkeypatch, tmp_path):
    reactome_path = _write_tsv(
        tmp_path / "reactome.tsv",
        pd.DataFrame(
            {
                "Pathway identifier": ["R-HSA-1"],
                "Pathway name": ["Test Pathway"],
                "Entity identifier": ["ENSG000001"],
            }
        ),
    )
    opentargets_path = _write_tsv(
        tmp_path / "opentargets.tsv",
        pd.DataFrame(
            {
                "diseaseId": ["EFO_0001"],
                "targetId": ["ENSG000001"],
                "overallScore": [0.8],
            }
        ),
    )
    gwas_path = _write_tsv(
        tmp_path / "gwas.tsv",
        pd.DataFrame(
            {
                "trait_id": ["EFO_0001"],
                "gene_id": ["ENSG000002"],
                "p_value": [1e-6],
            }
        ),
    )
    disgenet_path = _write_tsv(
        tmp_path / "disgenet.tsv",
        pd.DataFrame(
            {
                "diseaseId": ["EFO_0001"],
                "geneId": ["ENSG000003"],
                "score": [0.2],
            }
        ),
    )
    hpa_path = _write_tsv(
        tmp_path / "hpa.tsv",
        pd.DataFrame(
            {
                "Gene": ["ENSG000001"],
                "Tissue": ["Retina"],
                "TPM": [50.0],
            }
        ),
    )
    plae_path = _write_tsv(
        tmp_path / "plae.tsv",
        pd.DataFrame(
            {
                "gene_id": ["ENSG000001"],
                "tissue": ["Retina"],
                "log2_tpm": [5.5],
            }
        ),
    )
    rxnorm_path = _write_tsv(
        tmp_path / "rxnorm.tsv",
        pd.DataFrame(
            {
                "RXCUI": ["123"],
                "STR": ["TestDrug"],
                "synonym": ["Drug Synonym"],
            }
        ),
    )
    sider_path = _write_tsv(
        tmp_path / "sider.tsv",
        pd.DataFrame(
            {
                "drug_id": ["123"],
                "adverse_event": ["Headache"],
                "reports": [10],
            }
        ),
    )
    trials_path = _write_dataframe(
        tmp_path / "trials.csv",
        pd.DataFrame(
            {
                "nct_id": ["NCT0001"],
                "brief_title": ["A Trial"],
            }
        ),
    )

    omni_df = pd.DataFrame(
        {
            "source_genesymbol": ["GENE1"],
            "target_genesymbol": ["GENE2"],
            "is_stimulation": [True],
            "references": ["PMID:1"],
        }
    )
    string_df = pd.DataFrame(
        {
            "protein1": ["9606.GENE1"],
            "protein2": ["9606.GENE3"],
            "combined_score": [600],
        }
    )
    drugcentral_df = pd.DataFrame(
        {
            "drug_chembl_id": ["CHEMBL1"],
            "gene": ["GENE1"],
            "action_type": ["inhibitor"],
            "act_value": [20.0],
            "act_unit": ["nM"],
            "reference": ["PMID:2"],
        }
    )
    iuphar_df = pd.DataFrame(
        {
            "ligand_id": ["CHEMBL1"],
            "uniprot_id": ["GENE2"],
            "mechanism_of_action": ["agonist"],
        }
    )

    monkeypatch.setattr(build_graph, "_load_omnipath", lambda: omni_df)
    monkeypatch.setattr(build_graph, "_load_signor", lambda: pd.DataFrame())
    monkeypatch.setattr(build_graph, "_load_string", lambda: string_df)
    monkeypatch.setattr(build_graph, "_load_drugcentral", lambda: drugcentral_df)
    monkeypatch.setattr(build_graph, "_load_iuphar", lambda: iuphar_df)

    config = build_graph.BuildConfig(
        sqlite_path=None,
        reactome=reactome_path,
        opentargets=opentargets_path,
        gwas=gwas_path,
        disgenet=disgenet_path,
        hpa=hpa_path,
        plae=plae_path,
        rxnorm=rxnorm_path,
        sider=sider_path,
        clinicaltrials=trials_path,
    )
    tables = build_graph.assemble_tables(config)

    assert set(tables.keys()) >= {
        "protein_edge",
        "pathway_member",
        "disease_gene",
        "tissue_expr",
        "drug",
        "drug_target",
        "safety_ae",
        "trial",
    }
    assert len(tables["protein_edge"]) == 2
    assert len(tables["drug_target"]) == 2
    assert len(tables["disease_gene"]) == 3
    assert len(tables["drug"]) == 1
    assert len(tables["trial"]) == 1


def test_build_knowledge_graph_sqlite(monkeypatch, tmp_path):
    monkeypatch.setattr(build_graph, "assemble_tables", lambda config: {
        "protein_edge": pd.DataFrame(
            [
                {
                    "src_gene_id": "GENE1",
                    "dst_gene_id": "GENE2",
                    "relation": "interaction",
                    "sign": "+",
                    "direct": True,
                    "evidence": "{}",
                    "source": "Test",
                    "source_reference": "PMID:1",
                }
            ]
        ),
        "pathway_member": pd.DataFrame(),
        "disease_gene": pd.DataFrame(),
        "tissue_expr": pd.DataFrame(),
        "drug": pd.DataFrame(),
        "drug_target": pd.DataFrame(),
        "safety_ae": pd.DataFrame(),
        "trial": pd.DataFrame(),
    })
    db_path = tmp_path / "kg.sqlite"
    config = build_graph.BuildConfig(sqlite_path=db_path)
    build_graph.build_knowledge_graph(config)
    assert db_path.exists()
    with sqlite3.connect(db_path) as conn:
        rows = conn.execute("SELECT COUNT(*) FROM protein_edge").fetchone()[0]
    assert rows == 1
