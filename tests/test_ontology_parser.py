import pytest

pd = pytest.importorskip("pandas")

from etl.ontology_parser import (
    DEFAULT_OCULAR_ROOTS,
    annotate_phenotypes,
    build_manifest,
)


def test_annotate_phenotypes_assigns_scopes():
    phenotype_df = pd.DataFrame(
        [
            {"id": "HP:root_eye", "name": "Eye abnormality"},
            {"id": "HP:child_eye", "name": "Retinal issue"},
            {"id": "HP:other", "name": "Non ocular"},
        ]
    )

    ocular_descendants = {
        "HP:root_eye": {"HP:root_eye", "HP:child_eye"},
        "HP:root_other": {"HP:other"},
    }

    annotated, rows = annotate_phenotypes(
        phenotype_df,
        ocular_descendants,
        ocular_roots={"HP:root_eye": "eye", "HP:root_other": "other"},
    )

    child_row = annotated.loc[annotated["id"] == "HP:child_eye"].iloc[0]
    assert child_row["is_ocular"] is True
    assert child_row["ocular_scope"] == "eye"

    other_row = annotated.loc[annotated["id"] == "HP:other"].iloc[0]
    assert other_row["is_ocular"] is True
    assert other_row["ocular_scope"] == "other"

    assert {entry["root_scope"] for entry in rows} == {"eye", "other"}


def test_build_manifest_counts_records():
    disease_df = pd.DataFrame({"id": ["D1", "D2"]})
    phenotype_df = pd.DataFrame(
        {
            "id": ["HP:1", "HP:2"],
            "is_ocular": [True, False],
        }
    )
    bridge_df = pd.DataFrame({"disease_id": ["D1"], "phenotype_id": ["HP:1"]})
    ocular_descendants = {key: {key} for key in DEFAULT_OCULAR_ROOTS.keys()}

    manifest = build_manifest(
        "2024-01-01",
        disease_df,
        phenotype_df,
        bridge_df,
        ocular_descendants,
        source_files={"mondo": "mondo.obo"},
    )

    assert manifest["release"] == "2024-01-01"
    assert manifest["disease_records"] == 2
    assert manifest["phenotype_records"] == 2
    assert manifest["ocular_phenotype_records"] == 1
    assert manifest["disease_phenotype_records"] == 1
    assert manifest["source_files"]["mondo"] == "mondo.obo"
