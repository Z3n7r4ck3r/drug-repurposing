from __future__ import annotations

import pytest

pytest.importorskip("pandas")

from graph import schema


def test_schema_sql_contains_tables(tmp_path):
    ddl_path = tmp_path / "schema.sql"
    schema.write_schema_sql(ddl_path)
    content = ddl_path.read_text()
    assert "CREATE TABLE IF NOT EXISTS gene" in content
    assert "CREATE TABLE IF NOT EXISTS drug_target" in content


def test_schema_manifest_columns():
    manifest = schema.schema_manifest()
    assert "protein_edge" in manifest
    assert manifest["protein_edge"]["src_gene_id"] == "TEXT"
    assert manifest["trial"]["nct_id"] == "TEXT"
