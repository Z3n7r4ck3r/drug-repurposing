"""Loaders for genetics evidence (OpenTargets, GWAS, DisGeNET)."""
from __future__ import annotations

import json
from pathlib import Path
from typing import List, Mapping, Optional

import pandas as pd

DISEASE_GENE_COLUMNS = ["disease_id", "gene_id", "evidence_type", "score", "source", "evidence"]


def _read_table(path: Optional[str | Path]) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame()
    file_path = Path(path)
    if not file_path.exists():
        return pd.DataFrame()
    sep = "\t" if file_path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(file_path, sep=sep)


def load_opentargets(path: Optional[str | Path]) -> pd.DataFrame:
    df = _read_table(path)
    if df.empty:
        return pd.DataFrame(columns=DISEASE_GENE_COLUMNS)

    records: List[Mapping[str, object]] = []
    for record in df.to_dict(orient="records"):
        disease = record.get("diseaseId") or record.get("disease_id")
        gene = record.get("targetId") or record.get("gene_id")
        if not disease or not gene:
            continue
        score = record.get("overallScore") or record.get("score")
        datatype = record.get("datatypeId") or record.get("data_source") or "opentargets"
        evidence = {"data_type": datatype}
        association_score = record.get("association_score")
        if association_score is not None:
            evidence["association_score"] = association_score
        records.append(
            {
                "disease_id": str(disease),
                "gene_id": str(gene),
                "evidence_type": str(datatype),
                "score": float(score) if pd.notna(score) else None,
                "source": "OpenTargets",
                "evidence": json.dumps(evidence),
            }
        )

    result = pd.DataFrame(records, columns=DISEASE_GENE_COLUMNS)
    result.drop_duplicates(subset=["disease_id", "gene_id", "evidence_type"], inplace=True)
    return result


def load_gwas(path: Optional[str | Path]) -> pd.DataFrame:
    df = _read_table(path)
    if df.empty:
        return pd.DataFrame(columns=DISEASE_GENE_COLUMNS)

    records: List[Mapping[str, object]] = []
    for record in df.to_dict(orient="records"):
        disease = record.get("trait_id") or record.get("disease_id") or record.get("efo_id")
        gene = record.get("gene_id") or record.get("ensembl_id")
        pval = record.get("p_value") or record.get("pvalue")
        odds = record.get("odds_ratio") or record.get("or")
        if not disease or not gene:
            continue
        evidence = {}
        if pval is not None and pd.notna(pval):
            evidence["p_value"] = float(pval)
        if odds is not None and pd.notna(odds):
            evidence["odds_ratio"] = float(odds)
        records.append(
            {
                "disease_id": str(disease),
                "gene_id": str(gene),
                "evidence_type": "GWAS",
                "score": float(record.get("beta")) if record.get("beta") is not None else None,
                "source": "GWASCatalog",
                "evidence": json.dumps(evidence) if evidence else json.dumps({}),
            }
        )
    result = pd.DataFrame(records, columns=DISEASE_GENE_COLUMNS)
    result.drop_duplicates(subset=["disease_id", "gene_id", "evidence_type"], inplace=True)
    return result


def load_disgenet(path: Optional[str | Path]) -> pd.DataFrame:
    df = _read_table(path)
    if df.empty:
        return pd.DataFrame(columns=DISEASE_GENE_COLUMNS)

    records: List[Mapping[str, object]] = []
    for record in df.to_dict(orient="records"):
        disease = record.get("diseaseId") or record.get("disease_id")
        gene = record.get("geneId") or record.get("gene_id")
        if not disease or not gene:
            continue
        score = record.get("score") or record.get("DSI") or record.get("DPI")
        evidence = {"source": record.get("source")}
        records.append(
            {
                "disease_id": str(disease),
                "gene_id": str(gene),
                "evidence_type": "DisGeNET",
                "score": float(score) if score is not None and pd.notna(score) else None,
                "source": "DisGeNET",
                "evidence": json.dumps(evidence),
            }
        )
    result = pd.DataFrame(records, columns=DISEASE_GENE_COLUMNS)
    result.drop_duplicates(subset=["disease_id", "gene_id", "evidence_type"], inplace=True)
    return result


__all__ = ["load_opentargets", "load_gwas", "load_disgenet", "DISEASE_GENE_COLUMNS"]
