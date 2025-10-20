"""Build a normalized knowledge graph and load it into SQLite/PostgreSQL."""
from __future__ import annotations

import argparse
import json
import logging
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence

import pandas as pd

from etl.clinicaltrials_etl import load_clinicaltrials
from etl.expression_etl import load_hpa, load_plae
from etl.genetics_etl import load_disgenet, load_gwas, load_opentargets
from etl.reactome_etl import load_reactome
from etl.rxnorm_etl import load_rxnorm_drugs
from etl.sider_etl import load_sider
from graph import schema as schema_module

LOGGER = logging.getLogger(__name__)

EDGE_COLUMNS = [
    "src_gene_id",
    "dst_gene_id",
    "relation",
    "sign",
    "direct",
    "evidence",
    "source",
    "source_reference",
]

DRUG_TARGET_COLUMNS = [
    "drug_id",
    "target_id",
    "target_type",
    "action",
    "affinity",
    "affinity_unit",
    "moa_category",
    "source",
    "evidence",
]

PATHWAY_COLUMNS = ["pathway_id", "gene_id", "pathway_name", "source", "evidence"]
DISEASE_GENE_COLUMNS = ["disease_id", "gene_id", "evidence_type", "score", "source", "evidence"]
TISSUE_COLUMNS = ["gene_id", "tissue", "expression", "unit", "source", "evidence"]
DRUG_COLUMNS = ["drug_id", "preferred_name", "synonyms", "source"]
SAFETY_COLUMNS = [
    "drug_id",
    "adverse_event",
    "report_count",
    "proportional_reporting_ratio",
    "source",
    "evidence",
]
TRIAL_COLUMNS = [
    "nct_id",
    "title",
    "status",
    "phase",
    "conditions",
    "enrollment",
    "interventions",
    "last_updated",
    "source",
]


def _truthy(value) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    if isinstance(value, (int, float)):
        return value != 0
    text = str(value).strip().lower()
    return text in {"1", "true", "t", "yes", "y", "direct", "up"}


def _first(record: Mapping, keys: Iterable[str]) -> Optional[str]:
    for key in keys:
        if key not in record:
            continue
        value = record.get(key)
        if value is None:
            continue
        if isinstance(value, float) and pd.isna(value):
            continue
        text = str(value).strip()
        if text:
            return text
    return None


def _safe_fetch(fetcher, name: str) -> pd.DataFrame:
    try:
        return fetcher()
    except Exception as exc:  # pragma: no cover - optional connectors may fail
        LOGGER.warning("Skipping %s due to error: %s", name, exc)
        return pd.DataFrame()


def _load_omnipath() -> pd.DataFrame:
    try:
        from etl.omnipath_etl import get_omnipath_interactions  # type: ignore
    except Exception as exc:  # pragma: no cover - optional dependency
        LOGGER.warning("OmniPath connector unavailable: %s", exc)
        return pd.DataFrame()
    return _safe_fetch(get_omnipath_interactions, "OmniPath")


def _load_signor() -> pd.DataFrame:
    try:
        from etl.signor_etl import get_signor_interactions  # type: ignore
    except Exception as exc:  # pragma: no cover - optional dependency
        LOGGER.warning("SIGNOR connector unavailable: %s", exc)
        return pd.DataFrame()
    return _safe_fetch(get_signor_interactions, "SIGNOR")


def _load_string() -> pd.DataFrame:
    try:
        from etl.string_etl import get_string_interactions  # type: ignore
    except Exception as exc:  # pragma: no cover - optional dependency
        LOGGER.warning("STRING connector unavailable: %s", exc)
        return pd.DataFrame()
    return _safe_fetch(get_string_interactions, "STRING")


def _load_drugcentral() -> pd.DataFrame:
    try:
        from etl.drugcentral_etl import get_drugcentral_interactions  # type: ignore
    except Exception as exc:  # pragma: no cover - optional dependency
        LOGGER.warning("DrugCentral connector unavailable: %s", exc)
        return pd.DataFrame()
    return _safe_fetch(get_drugcentral_interactions, "DrugCentral")


def _load_iuphar() -> pd.DataFrame:
    try:
        from etl.iuphar_etl import get_iuphar_interactions  # type: ignore
    except Exception as exc:  # pragma: no cover - optional dependency
        LOGGER.warning("IUPHAR connector unavailable: %s", exc)
        return pd.DataFrame()
    return _safe_fetch(get_iuphar_interactions, "IUPHAR")


def normalize_signed_edges(df: pd.DataFrame, source_name: str) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=EDGE_COLUMNS)

    records: List[Mapping[str, object]] = []
    for record in df.to_dict(orient="records"):
        src = _first(record, ["source_genesymbol", "source", "ENTITYA", "ENTITY_A", "ENTITYA_NAME"])
        dst = _first(record, ["target_genesymbol", "target", "ENTITYB", "ENTITY_B", "ENTITYB_NAME"])
        if not src or not dst:
            continue

        relation = (
            _first(record, ["consensus_direction", "interaction_type", "type", "mechanism", "MECHANISM", "EFFECT"])
            or "interaction"
        )
        relation = relation.replace(" ", "_")

        sign = None
        if any(
            _truthy(record.get(key))
            for key in ("is_stimulation", "consensus_stimulation", "stimulation", "UP_REGULATION", "up-regulates")
        ):
            sign = "+"
        elif any(
            _truthy(record.get(key))
            for key in ("is_inhibition", "consensus_inhibition", "inhibition", "DOWN_REGULATION", "down-regulates")
        ):
            sign = "-"
        effect = (_first(record, ["EFFECT"]) or "").lower()
        if sign is None:
            if "activ" in effect:
                sign = "+"
            elif "inhib" in effect or "down" in effect:
                sign = "-"
        if sign is None:
            continue

        direct = any(
            _truthy(record.get(key))
            for key in ("is_direct", "direct", "DIRECT", "is_directed")
        )
        evidence = _first(record, ["references", "curation_effort", "pmid", "REFERENCE", "PMID"])

        records.append(
            {
                "src_gene_id": src,
                "dst_gene_id": dst,
                "relation": relation,
                "sign": sign,
                "direct": bool(direct or source_name == "OmniPath"),
                "evidence": json.dumps({"raw": evidence}) if evidence else json.dumps({}),
                "source": source_name,
                "source_reference": evidence or "",
            }
        )

    edges = pd.DataFrame(records, columns=EDGE_COLUMNS)
    edges.drop_duplicates(subset=["src_gene_id", "dst_gene_id", "relation", "sign", "source"], inplace=True)
    return edges


def normalize_string(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=EDGE_COLUMNS)

    records: List[Mapping[str, object]] = []
    for record in df.to_dict(orient="records"):
        src = record.get("protein1")
        dst = record.get("protein2")
        if not src or not dst:
            continue
        if isinstance(src, str) and "." in src:
            src = src.split(".")[-1]
        if isinstance(dst, str) and "." in dst:
            dst = dst.split(".")[-1]
        score = record.get("combined_score")
        evidence = json.dumps({"combined_score": float(score)}) if score is not None else json.dumps({})
        records.append(
            {
                "src_gene_id": str(src),
                "dst_gene_id": str(dst),
                "relation": "physical_interaction",
                "sign": None,
                "direct": False,
                "evidence": evidence,
                "source": "STRING",
                "source_reference": "STRING",
            }
        )

    edges = pd.DataFrame(records, columns=EDGE_COLUMNS)
    edges.drop_duplicates(subset=["src_gene_id", "dst_gene_id", "source"], inplace=True)
    return edges


def normalize_drug_targets(df: pd.DataFrame, source_name: str) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=DRUG_TARGET_COLUMNS)

    records: List[Mapping[str, object]] = []
    for record in df.to_dict(orient="records"):
        drug_id = _first(
            record,
            [
                "drug_chembl_id",
                "molecule_chembl_id",
                "drugcentral_id",
                "drugbank_id",
                "ligand_id",
                "drug_id",
                "DRUG_ID",
            ],
        )
        target_id = _first(
            record,
            [
                "gene",
                "target_gene_symbol",
                "accession",
                "uniprot_id",
                "swissprot",
                "target_id",
                "target",
            ],
        )
        if not drug_id or not target_id:
            continue

        action = _first(
            record,
            ["action_type", "action", "act_comment", "activity_comment", "mechanism_of_action", "mode_of_action"],
        )
        affinity_value = _first(record, ["act_value", "standard_value", "affinity", "pchembl_value", "activity_value"])
        affinity_unit = _first(record, ["act_unit", "standard_units", "affinity_unit"])
        target_type = _first(
            record,
            ["target_class", "target_type", "target_pref_name", "target_organism"],
        )
        moa = _first(record, ["moa", "mechanism_comment", "mechanism_of_action"])
        references = _first(record, ["reference", "pmid", "pubmed_id", "source_reference"])

        evidence = {}
        if references:
            evidence["references"] = references
        if record.get("relation"):
            evidence["relation"] = record["relation"]

        affinity_numeric = None
        if affinity_value is not None:
            try:
                affinity_numeric = float(affinity_value)
            except ValueError:
                affinity_numeric = None

        records.append(
            {
                "drug_id": drug_id,
                "target_id": target_id,
                "target_type": target_type or "",
                "action": action or "",
                "affinity": affinity_numeric,
                "affinity_unit": affinity_unit or "",
                "moa_category": moa or "",
                "source": source_name,
                "evidence": json.dumps(evidence) if evidence else json.dumps({}),
            }
        )

    targets = pd.DataFrame(records, columns=DRUG_TARGET_COLUMNS)
    targets.drop_duplicates(subset=["drug_id", "target_id", "source"], inplace=True)
    return targets


@dataclass
class BuildConfig:
    sqlite_path: Optional[Path] = None
    postgres_dsn: Optional[str] = None
    reactome: Optional[Path] = None
    opentargets: Optional[Path] = None
    gwas: Optional[Path] = None
    disgenet: Optional[Path] = None
    hpa: Optional[Path] = None
    plae: Optional[Path] = None
    rxnorm: Optional[Path] = None
    sider: Optional[Path] = None
    clinicaltrials: Optional[Path] = None


def assemble_tables(config: BuildConfig) -> Dict[str, pd.DataFrame]:
    LOGGER.info("Loading causal interaction datasets")
    omnipath_edges = normalize_signed_edges(_load_omnipath(), "OmniPath")
    signor_edges = normalize_signed_edges(_load_signor(), "SIGNOR")
    string_edges = normalize_string(_load_string())
    protein_edges = pd.concat([omnipath_edges, signor_edges, string_edges], ignore_index=True)
    protein_edges.drop_duplicates(
        subset=["src_gene_id", "dst_gene_id", "relation", "sign", "source"], inplace=True
    )
    if protein_edges.empty:
        protein_edges = pd.DataFrame(columns=EDGE_COLUMNS)

    LOGGER.info("Loading pathway memberships")
    reactome = load_reactome(config.reactome)
    if reactome.empty:
        reactome = pd.DataFrame(columns=PATHWAY_COLUMNS)
    else:
        reactome = reactome.assign(source="Reactome")

    LOGGER.info("Loading genetics evidence")
    genetics_frames = [
        load_opentargets(config.opentargets),
        load_gwas(config.gwas),
        load_disgenet(config.disgenet),
    ]
    disease_gene = pd.concat(genetics_frames, ignore_index=True)
    if disease_gene.empty:
        disease_gene = pd.DataFrame(columns=DISEASE_GENE_COLUMNS)

    LOGGER.info("Loading expression datasets")
    hpa = load_hpa(config.hpa)
    plae = load_plae(config.plae)
    tissue_expr = pd.concat([hpa, plae], ignore_index=True)
    if tissue_expr.empty:
        tissue_expr = pd.DataFrame(columns=TISSUE_COLUMNS)

    LOGGER.info("Loading drug catalogues")
    rxnorm = load_rxnorm_drugs(config.rxnorm)
    drugcentral_targets = normalize_drug_targets(_load_drugcentral(), "DrugCentral")
    iuphar_targets = normalize_drug_targets(_load_iuphar(), "IUPHAR")
    drug_target = pd.concat([drugcentral_targets, iuphar_targets], ignore_index=True)
    if drug_target.empty:
        drug_target = pd.DataFrame(columns=DRUG_TARGET_COLUMNS)

    LOGGER.info("Loading safety and trials data")
    safety = load_sider(config.sider)
    trials = load_clinicaltrials(config.clinicaltrials)

    tables: Dict[str, pd.DataFrame] = {
        "protein_edge": protein_edges,
        "pathway_member": reactome,
        "disease_gene": disease_gene,
        "tissue_expr": tissue_expr,
        "drug": rxnorm if not rxnorm.empty else pd.DataFrame(columns=DRUG_COLUMNS),
        "drug_target": drug_target,
        "safety_ae": safety if not safety.empty else pd.DataFrame(columns=SAFETY_COLUMNS),
        "trial": trials if not trials.empty else pd.DataFrame(columns=TRIAL_COLUMNS),
    }
    return tables


SQLITE_DDL = {
    "protein_edge": """
        CREATE TABLE IF NOT EXISTS protein_edge (
            edge_id INTEGER PRIMARY KEY AUTOINCREMENT,
            src_gene_id TEXT NOT NULL,
            dst_gene_id TEXT NOT NULL,
            relation TEXT NOT NULL,
            sign TEXT,
            direct INTEGER DEFAULT 0,
            evidence TEXT,
            source TEXT NOT NULL,
            source_reference TEXT
        )
    """,
    "pathway_member": """
        CREATE TABLE IF NOT EXISTS pathway_member (
            pathway_id TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            pathway_name TEXT,
            source TEXT NOT NULL,
            evidence TEXT
        )
    """,
    "disease_gene": """
        CREATE TABLE IF NOT EXISTS disease_gene (
            disease_id TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            evidence_type TEXT NOT NULL,
            score REAL,
            source TEXT NOT NULL,
            evidence TEXT
        )
    """,
    "tissue_expr": """
        CREATE TABLE IF NOT EXISTS tissue_expr (
            record_id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_id TEXT NOT NULL,
            tissue TEXT NOT NULL,
            expression REAL,
            unit TEXT,
            source TEXT NOT NULL,
            evidence TEXT
        )
    """,
    "drug": """
        CREATE TABLE IF NOT EXISTS drug (
            drug_id TEXT PRIMARY KEY,
            preferred_name TEXT NOT NULL,
            synonyms TEXT,
            source TEXT NOT NULL
        )
    """,
    "drug_target": """
        CREATE TABLE IF NOT EXISTS drug_target (
            drug_id TEXT NOT NULL,
            target_id TEXT NOT NULL,
            target_type TEXT,
            action TEXT,
            affinity REAL,
            affinity_unit TEXT,
            moa_category TEXT,
            source TEXT NOT NULL,
            evidence TEXT
        )
    """,
    "safety_ae": """
        CREATE TABLE IF NOT EXISTS safety_ae (
            drug_id TEXT NOT NULL,
            adverse_event TEXT NOT NULL,
            report_count INTEGER,
            proportional_reporting_ratio REAL,
            source TEXT NOT NULL,
            evidence TEXT
        )
    """,
    "trial": """
        CREATE TABLE IF NOT EXISTS trial (
            nct_id TEXT PRIMARY KEY,
            title TEXT NOT NULL,
            status TEXT,
            phase TEXT,
            conditions TEXT,
            enrollment INTEGER,
            interventions TEXT,
            last_updated TEXT,
            source TEXT NOT NULL
        )
    """,
}


def _initialise_sqlite(conn: sqlite3.Connection) -> None:
    for statement in SQLITE_DDL.values():
        conn.execute(statement)
    conn.commit()


def _write_sqlite(conn: sqlite3.Connection, tables: Dict[str, pd.DataFrame]) -> None:
    for table, df in tables.items():
        if table not in SQLITE_DDL:
            LOGGER.debug("Skipping unknown table %s for SQLite load", table)
            continue
        if df.empty:
            LOGGER.info("No rows to write for %s", table)
            continue
        LOGGER.info("Writing %s rows to %s", len(df), table)
        df.to_sql(table, conn, if_exists="append", index=False)
    conn.commit()


def _write_postgres(dsn: str, tables: Dict[str, pd.DataFrame]) -> None:
    try:
        import psycopg
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("psycopg is required for Postgres loading") from exc

    LOGGER.info("Applying schema to PostgreSQL instance")
    with psycopg.connect(dsn, autocommit=True) as conn:  # pragma: no cover - requires postgres
        with conn.cursor() as cursor:
            schema_module.apply_schema(cursor)

    LOGGER.info("Loading tables into PostgreSQL")
    with psycopg.connect(dsn) as conn:  # pragma: no cover - requires postgres
        for table, df in tables.items():
            if df.empty:
                continue
            columns = list(df.columns)
            placeholders = ",".join(["%s"] * len(columns))
            insert_sql = f"INSERT INTO {table} ({','.join(columns)}) VALUES ({placeholders})"
            records = [tuple(row) for row in df.itertuples(index=False, name=None)]
            with conn.cursor() as cursor:
                cursor.execute(f"TRUNCATE {table}")
                for chunk_start in range(0, len(records), 1000):
                    chunk = records[chunk_start : chunk_start + 1000]
                    cursor.executemany(insert_sql, chunk)
        conn.commit()


def build_knowledge_graph(config: BuildConfig) -> Dict[str, pd.DataFrame]:
    tables = assemble_tables(config)

    if config.sqlite_path:
        config.sqlite_path.parent.mkdir(parents=True, exist_ok=True)
        LOGGER.info("Writing SQLite database to %s", config.sqlite_path)
        with sqlite3.connect(config.sqlite_path) as conn:
            _initialise_sqlite(conn)
            _write_sqlite(conn, tables)

    if config.postgres_dsn:
        LOGGER.info("Loading tables into PostgreSQL: %s", config.postgres_dsn)
        _write_postgres(config.postgres_dsn, tables)

    return tables


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sqlite-path", type=Path, help="Path to the SQLite database to materialise")
    parser.add_argument("--postgres-dsn", type=str, help="Optional psycopg DSN for a PostgreSQL target")
    parser.add_argument("--reactome", type=Path, help="Path to Reactome membership file")
    parser.add_argument("--opentargets", type=Path, help="Path to OpenTargets association export")
    parser.add_argument("--gwas", type=Path, help="Path to GWAS summary association export")
    parser.add_argument("--disgenet", type=Path, help="Path to DisGeNET association export")
    parser.add_argument("--hpa", type=Path, help="Path to HPA expression data")
    parser.add_argument("--plae", type=Path, help="Path to PLAE expression data")
    parser.add_argument("--rxnorm", type=Path, help="Path to RxNorm concepts export")
    parser.add_argument("--sider", type=Path, help="Path to SIDER adverse event export")
    parser.add_argument("--clinicaltrials", type=Path, help="Path to ClinicalTrials.gov export")
    parser.add_argument("--log-level", default="INFO")
    return parser


def main(argv: Sequence[str] | None = None) -> None:  # pragma: no cover - CLI wrapper
    parser = _build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO))

    config = BuildConfig(
        sqlite_path=getattr(args, "sqlite_path", None),
        postgres_dsn=getattr(args, "postgres_dsn", None),
        reactome=getattr(args, "reactome", None),
        opentargets=getattr(args, "opentargets", None),
        gwas=getattr(args, "gwas", None),
        disgenet=getattr(args, "disgenet", None),
        hpa=getattr(args, "hpa", None),
        plae=getattr(args, "plae", None),
        rxnorm=getattr(args, "rxnorm", None),
        sider=getattr(args, "sider", None),
        clinicaltrials=getattr(args, "clinicaltrials", None),
    )
    build_knowledge_graph(config)


if __name__ == "__main__":
    main()
