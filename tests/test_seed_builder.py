import pytest

pd = pytest.importorskip("pandas")

from signatures.seed_builder import SeedEvidence, assemble_seed_table


def test_assemble_seed_table_merges_sources():
    ot = pd.DataFrame(
        {
            "disease_id": ["D1"],
            "gene_symbol": ["GENE1"],
            "score": [0.5],
            "source": ["OpenTargets"],
            "evidence_type": ["genetic"],
        }
    )
    gwas = pd.DataFrame(
        {
            "disease_id": ["D1"],
            "gene_symbol": ["GENE2"],
            "score": [1.0],
            "source": ["GWAS"],
            "evidence_type": ["genetic"],
        }
    )
    curated = [SeedEvidence("D1", "GENE3", 0.8, "Literature", "curated", "PMID:1")]

    seeds = assemble_seed_table(ot, gwas, curated)
    assert len(seeds) == 3
    assert set(seeds["gene_symbol"]) == {"GENE1", "GENE2", "GENE3"}
