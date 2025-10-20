"""Utilities to build disease gene seed sets from heterogeneous evidence sources."""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Iterable, List, Optional

import numpy as np
import pandas as pd


@dataclass
class SeedEvidence:
    """Representation of a single disease→gene association."""

    disease_id: str
    gene_symbol: str
    score: float
    source: str
    evidence_type: str
    details: Optional[str] = None

    def as_dict(self) -> dict:
        payload = {
            "disease_id": self.disease_id,
            "gene_symbol": self.gene_symbol,
            "score": float(self.score),
            "source": self.source,
            "evidence_type": self.evidence_type,
        }
        if self.details:
            payload["details"] = self.details
        return payload


def _safe_fetch(callable_obj, name: str) -> pd.DataFrame:
    try:
        return callable_obj()
    except Exception as exc:  # pragma: no cover - defensive around optional deps
        logging.warning("Skipping %s seeds due to error: %s", name, exc)
        return pd.DataFrame()


def load_opentargets_seeds(min_score: float = 0.2) -> pd.DataFrame:
    """Fetch and normalise OpenTargets loci-to-gene scores."""

    try:
        from etl.opentargets_etl import get_opentargets_associations  # type: ignore
    except Exception as exc:  # pragma: no cover - optional dependency
        logging.warning("OpenTargets connector unavailable: %s", exc)
        return pd.DataFrame(columns=["disease_id", "gene_symbol", "score", "source", "evidence_type", "details"])

    raw = _safe_fetch(get_opentargets_associations, "OpenTargets")
    if raw.empty:
        return pd.DataFrame(columns=["disease_id", "gene_symbol", "score", "source", "evidence_type", "details"])

    disease_col = next((col for col in ["study_locus_id", "disease_id", "trait_id"] if col in raw.columns), None)
    gene_col = next((col for col in ["gene_id", "gene_symbol", "target_gene_id"] if col in raw.columns), None)
    score_col = next((col for col in ["overall_score", "score", "posterior_prob"] if col in raw.columns), None)

    if not disease_col or not gene_col or not score_col:
        logging.warning("OpenTargets dataframe missing key columns; got columns: %s", list(raw.columns))
        return pd.DataFrame(columns=["disease_id", "gene_symbol", "score", "source", "evidence_type", "details"])

    df = raw[[disease_col, gene_col, score_col]].rename(
        columns={disease_col: "disease_id", gene_col: "gene_symbol", score_col: "score"}
    )
    df = df.dropna(subset=["disease_id", "gene_symbol", "score"])
    df = df[df["score"] >= float(min_score)]
    df["source"] = "OpenTargets"
    df["evidence_type"] = "genetic"
    return df


def load_gwas_seeds(efo_ids: Iterable[str], max_pvalue: float = 5e-8) -> pd.DataFrame:
    """Query the GWAS catalog associations for selected traits."""

    try:
        from etl.gwas_catalog_etl import get_gwas_catalog_associations_for_disease  # type: ignore
    except Exception as exc:  # pragma: no cover
        logging.warning("GWAS connector unavailable: %s", exc)
        return pd.DataFrame(columns=["disease_id", "gene_symbol", "score", "source", "evidence_type", "details"])

    frames: List[pd.DataFrame] = []
    for efo_id in efo_ids:
        associations = _safe_fetch(lambda eid=efo_id: get_gwas_catalog_associations_for_disease(eid), f"GWAS {efo_id}")
        if associations.empty:
            continue
        gene_col = next((col for col in ["reportedGene", "gene", "mappedGenes"] if col in associations.columns), None)
        pvalue_col = next((col for col in ["pValueMantissa", "pvalue", "pValue"] if col in associations.columns), None)
        exponent_col = "pValueExponent" if "pValueExponent" in associations.columns else None

        if not gene_col or not pvalue_col:
            logging.warning("GWAS dataframe for %s missing expected columns", efo_id)
            continue

        df = associations[[gene_col]].copy()
        df.rename(columns={gene_col: "gene_symbol"}, inplace=True)
        df["gene_symbol"] = df["gene_symbol"].astype(str)

        if exponent_col:
            df["pvalue"] = (associations[pvalue_col].astype(float) * (10 ** associations[exponent_col].astype(float))).astype(float)
        else:
            df["pvalue"] = associations[pvalue_col].astype(float)

        df = df[df["pvalue"] <= max_pvalue]
        if df.empty:
            continue
        df["disease_id"] = efo_id
        df["score"] = -df["pvalue"].apply(lambda x: pd.NA if x <= 0 else np.log10(x))
        df["score"].fillna(-8, inplace=True)
        df["source"] = "GWAS"
        df["evidence_type"] = "genetic"
        frames.append(df[["disease_id", "gene_symbol", "score", "source", "evidence_type"]])

    if not frames:
        return pd.DataFrame(columns=["disease_id", "gene_symbol", "score", "source", "evidence_type"])

    return pd.concat(frames, ignore_index=True)


def assemble_seed_table(
    opentargets: Optional[pd.DataFrame] = None,
    gwas: Optional[pd.DataFrame] = None,
    curated: Optional[Iterable[SeedEvidence]] = None,
) -> pd.DataFrame:
    """Merge the available seed evidence into a harmonised table."""

    frames: List[pd.DataFrame] = []
    if opentargets is not None and not opentargets.empty:
        frames.append(opentargets.copy())
    if gwas is not None and not gwas.empty:
        frames.append(gwas.copy())
    if curated:
        frames.append(pd.DataFrame([item.as_dict() for item in curated]))

    if not frames:
        return pd.DataFrame(columns=["disease_id", "gene_symbol", "score", "source", "evidence_type", "details"])

    merged = pd.concat(frames, ignore_index=True, sort=False)
    if "details" not in merged.columns:
        merged["details"] = ""
    merged["details"].fillna("", inplace=True)
    merged["score"] = pd.to_numeric(merged["score"], errors="coerce")
    merged = merged.dropna(subset=["disease_id", "gene_symbol", "score"])
    merged["disease_id"] = merged["disease_id"].astype(str)
    merged["gene_symbol"] = merged["gene_symbol"].astype(str)
    merged = merged.sort_values(["disease_id", "score"], ascending=[True, False])
    return merged.reset_index(drop=True)


def summarise_by_disease(seeds: pd.DataFrame) -> pd.DataFrame:
    """Return aggregate statistics (gene counts, mean score) per disease."""

    if seeds.empty:
        return pd.DataFrame(columns=["disease_id", "n_genes", "mean_score"])
    return (
        seeds.groupby("disease_id")
        .agg(n_genes=("gene_symbol", "nunique"), mean_score=("score", "mean"))
        .reset_index()
    )


def save_seed_table(seeds: pd.DataFrame, path: str) -> None:
    if seeds.empty:
        logging.warning("Seed table is empty; writing placeholder file to %s", path)
    seeds.to_csv(path, index=False)
    logging.info("Seed table with %s rows written to %s", len(seeds), path)
