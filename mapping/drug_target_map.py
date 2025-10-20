"""Utilities to harmonise drug target annotations across sources."""
from __future__ import annotations

import logging
import sqlite3
from pathlib import Path
from typing import Optional

import pandas as pd


def load_graph_targets(db_path: Path) -> pd.DataFrame:
    """Load drug-target relationships from the knowledge graph database."""

    if not db_path.exists():
        raise FileNotFoundError(f"Knowledge graph database not found: {db_path}")

    with sqlite3.connect(db_path) as conn:
        df = pd.read_sql_query("SELECT * FROM drug_target", conn)
    if df.empty:
        logging.warning("No drug_target records found in %s", db_path)
    return df


def attach_rxnorm_synonyms(targets: pd.DataFrame, rxnorm_table: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Attach RxNorm names to the provided drug target table."""

    if targets.empty:
        return targets.copy()

    if rxnorm_table is None:
        try:
            rxnorm_table = pd.read_csv("rxnorm_processed_data.csv")
        except FileNotFoundError:
            logging.warning("rxnorm_processed_data.csv not found; skipping synonym expansion")
            return targets.copy()

    rxnorm = rxnorm_table.rename(columns={"rxcui": "rxnorm_id", "str": "drug_name"})
    merged = targets.merge(rxnorm[["rxnorm_id", "drug_name"]], left_on="drug_id", right_on="rxnorm_id", how="left")
    merged["preferred_name"] = merged["drug_name"].fillna(merged["drug_id"])
    return merged


def summarise_by_drug(targets: pd.DataFrame) -> pd.DataFrame:
    if targets.empty:
        return pd.DataFrame(columns=["drug_id", "n_targets", "sources"])
    return (
        targets.groupby("drug_id")
        .agg(n_targets=("target_id", "nunique"), sources=("source", lambda x: "|".join(sorted(set(x)))))
        .reset_index()
    )


def export_mapping(targets: pd.DataFrame, path: Path) -> None:
    if targets.empty:
        logging.warning("Target mapping empty; writing placeholder file to %s", path)
    targets.to_csv(path, index=False)
    logging.info("Wrote %s drug-target associations to %s", len(targets), path)
