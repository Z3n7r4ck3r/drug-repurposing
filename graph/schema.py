"""PostgreSQL schema definition and utilities for the knowledge graph."""
from __future__ import annotations

import argparse
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class Column:
    name: str
    type: str
    constraints: Sequence[str] = ()

    def render(self) -> str:
        parts = [self.name, self.type]
        parts.extend(self.constraints)
        return " ".join(parts)


@dataclass(frozen=True)
class Table:
    name: str
    columns: Sequence[Column]
    indexes: Sequence[str] = ()

    def create_statement(self) -> str:
        column_sql = ",\n        ".join(col.render() for col in self.columns)
        return f"CREATE TABLE IF NOT EXISTS {self.name} (\n        {column_sql}\n);"


SCHEMA: Sequence[Table] = (
    Table(
        name="gene",
        columns=(
            Column("gene_id", "TEXT", ("PRIMARY KEY",)),
            Column("symbol", "TEXT", ("NOT NULL",)),
            Column("ensembl_id", "TEXT"),
            Column("uniprot_id", "TEXT"),
            Column("species", "TEXT", ("DEFAULT 'Homo sapiens'",)),
        ),
    ),
    Table(
        name="protein_edge",
        columns=(
            Column("edge_id", "BIGSERIAL", ("PRIMARY KEY",)),
            Column("src_gene_id", "TEXT", ("NOT NULL",)),
            Column("dst_gene_id", "TEXT", ("NOT NULL",)),
            Column("relation", "TEXT", ("NOT NULL",)),
            Column("sign", "TEXT"),
            Column("direct", "BOOLEAN", ("DEFAULT FALSE",)),
            Column("evidence", "TEXT"),
            Column("source", "TEXT", ("NOT NULL",)),
            Column("source_reference", "TEXT"),
        ),
        indexes=(
            "CREATE INDEX IF NOT EXISTS ix_protein_edge_src ON protein_edge(src_gene_id);",
            "CREATE INDEX IF NOT EXISTS ix_protein_edge_dst ON protein_edge(dst_gene_id);",
            "CREATE INDEX IF NOT EXISTS ix_protein_edge_source ON protein_edge(source);",
        ),
    ),
    Table(
        name="pathway_member",
        columns=(
            Column("pathway_id", "TEXT", ("NOT NULL",)),
            Column("gene_id", "TEXT", ("NOT NULL",)),
            Column("pathway_name", "TEXT"),
            Column("source", "TEXT", ("NOT NULL",)),
            Column("evidence", "TEXT"),
        ),
        indexes=(
            "CREATE INDEX IF NOT EXISTS ix_pathway_member_pathway ON pathway_member(pathway_id);",
            "CREATE INDEX IF NOT EXISTS ix_pathway_member_gene ON pathway_member(gene_id);",
        ),
    ),
    Table(
        name="disease_gene",
        columns=(
            Column("disease_id", "TEXT", ("NOT NULL",)),
            Column("gene_id", "TEXT", ("NOT NULL",)),
            Column("evidence_type", "TEXT", ("NOT NULL",)),
            Column("score", "DOUBLE PRECISION"),
            Column("source", "TEXT", ("NOT NULL",)),
            Column("evidence", "JSONB"),
        ),
        indexes=(
            "CREATE INDEX IF NOT EXISTS ix_disease_gene_disease ON disease_gene(disease_id);",
            "CREATE INDEX IF NOT EXISTS ix_disease_gene_gene ON disease_gene(gene_id);",
        ),
    ),
    Table(
        name="tissue_expr",
        columns=(
            Column("record_id", "BIGSERIAL", ("PRIMARY KEY",)),
            Column("gene_id", "TEXT", ("NOT NULL",)),
            Column("tissue", "TEXT", ("NOT NULL",)),
            Column("expression", "DOUBLE PRECISION"),
            Column("unit", "TEXT", ("DEFAULT 'TPM'",)),
            Column("source", "TEXT", ("NOT NULL",)),
        ),
        indexes=(
            "CREATE INDEX IF NOT EXISTS ix_tissue_expr_gene ON tissue_expr(gene_id);",
            "CREATE INDEX IF NOT EXISTS ix_tissue_expr_tissue ON tissue_expr(tissue);",
        ),
    ),
    Table(
        name="drug",
        columns=(
            Column("drug_id", "TEXT", ("PRIMARY KEY",)),
            Column("preferred_name", "TEXT", ("NOT NULL",)),
            Column("synonyms", "TEXT"),
            Column("smiles", "TEXT"),
            Column("inchikey", "TEXT"),
            Column("approval_status", "TEXT"),
            Column("source", "TEXT", ("NOT NULL",)),
        ),
    ),
    Table(
        name="drug_target",
        columns=(
            Column("drug_id", "TEXT", ("NOT NULL",)),
            Column("target_id", "TEXT", ("NOT NULL",)),
            Column("target_type", "TEXT"),
            Column("action", "TEXT"),
            Column("affinity", "DOUBLE PRECISION"),
            Column("affinity_unit", "TEXT"),
            Column("moa_category", "TEXT"),
            Column("source", "TEXT", ("NOT NULL",)),
            Column("evidence", "JSONB"),
        ),
        indexes=(
            "CREATE INDEX IF NOT EXISTS ix_drug_target_drug ON drug_target(drug_id);",
            "CREATE INDEX IF NOT EXISTS ix_drug_target_target ON drug_target(target_id);",
        ),
    ),
    Table(
        name="safety_ae",
        columns=(
            Column("drug_id", "TEXT", ("NOT NULL",)),
            Column("adverse_event", "TEXT", ("NOT NULL",)),
            Column("report_count", "INTEGER"),
            Column("proportional_reporting_ratio", "DOUBLE PRECISION"),
            Column("source", "TEXT", ("NOT NULL",)),
            Column("evidence", "JSONB"),
        ),
    ),
    Table(
        name="label_section",
        columns=(
            Column("drug_id", "TEXT", ("NOT NULL",)),
            Column("section", "TEXT", ("NOT NULL",)),
            Column("text", "TEXT", ("NOT NULL",)),
            Column("source", "TEXT", ("NOT NULL",)),
            Column("extracted_at", "TIMESTAMPTZ", ("DEFAULT NOW()",)),
        ),
    ),
    Table(
        name="trial",
        columns=(
            Column("nct_id", "TEXT", ("PRIMARY KEY",)),
            Column("title", "TEXT", ("NOT NULL",)),
            Column("status", "TEXT"),
            Column("phase", "TEXT"),
            Column("conditions", "TEXT"),
            Column("enrollment", "INTEGER"),
            Column("interventions", "TEXT"),
            Column("last_updated", "TIMESTAMPTZ"),
            Column("source", "TEXT", ("NOT NULL",)),
        ),
        indexes=(
            "CREATE INDEX IF NOT EXISTS ix_trial_status ON trial(status);",
            "CREATE INDEX IF NOT EXISTS ix_trial_phase ON trial(phase);",
        ),
    ),
)


def render_schema_sql() -> str:
    statements: List[str] = []
    for table in SCHEMA:
        statements.append(table.create_statement())
        statements.extend(table.indexes)
    return "\n\n".join(statements)


def apply_schema(cursor) -> None:
    """Apply the schema using a DB-API cursor.

    The function does not manage transactions; callers are responsible for
    committing. It works with psycopg (PostgreSQL) as well as sqlite3 for tests.
    """

    statements = render_schema_sql().split(";\n")
    for statement in statements:
        stmt = statement.strip()
        if not stmt:
            continue
        cursor.execute(stmt + ";")


def schema_manifest() -> Dict[str, Dict[str, str]]:
    return {
        table.name: {column.name: column.type for column in table.columns}
        for table in SCHEMA
    }


def write_schema_sql(path: Path) -> None:
    sql = render_schema_sql()
    path.write_text(sql + "\n", encoding="utf-8")
    LOGGER.info("Wrote schema DDL to %s", path)


def write_manifest(path: Path) -> None:
    manifest = schema_manifest()
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    LOGGER.info("Wrote schema manifest to %s", path)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Manage the knowledge graph schema")
    parser.add_argument(
        "--ddl-path",
        type=Path,
        default=Path("graph/schema.sql"),
        help="Where to write the schema SQL",
    )
    parser.add_argument(
        "--manifest-path",
        type=Path,
        default=Path("graph/schema_manifest.json"),
        help="Where to write the schema manifest",
    )
    parser.add_argument(
        "--apply-dsn",
        type=str,
        default=None,
        help="Optional psycopg-compatible DSN for applying the schema",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO)

    write_schema_sql(args.ddl_path)
    write_manifest(args.manifest_path)

    if args.apply_dsn:
        try:
            import psycopg
        except Exception as exc:  # pragma: no cover - optional dependency
            raise SystemExit("psycopg is required to apply schema: %s" % exc)

        LOGGER.info("Applying schema to %s", args.apply_dsn)
        with psycopg.connect(args.apply_dsn, autocommit=True) as conn:  # pragma: no cover - requires postgres
            with conn.cursor() as cursor:
                apply_schema(cursor)
        LOGGER.info("Schema applied successfully")


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    main()
