"""Loaders for tissue expression datasets (HPA, PLAE)."""
from __future__ import annotations

import json
from pathlib import Path
from typing import List, Mapping, Optional

import pandas as pd

TISSUE_COLUMNS = ["gene_id", "tissue", "expression", "unit", "source", "evidence"]


def _read_expression(path: Optional[str | Path]) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame()
    file_path = Path(path)
    if not file_path.exists():
        return pd.DataFrame()
    if file_path.suffix.lower() in {".parquet", ".pq"}:
        return pd.read_parquet(file_path)
    sep = "\t" if file_path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(file_path, sep=sep)


def load_hpa(path: Optional[str | Path]) -> pd.DataFrame:
    df = _read_expression(path)
    if df.empty:
        return pd.DataFrame(columns=TISSUE_COLUMNS)

    records: List[Mapping[str, object]] = []
    for record in df.to_dict(orient="records"):
        gene = record.get("Gene") or record.get("gene_id") or record.get("Ensembl")
        tissue = record.get("Tissue") or record.get("tissue")
        level = record.get("TPM") or record.get("expression")
        if not gene or not tissue:
            continue
        unit = "TPM"
        if "unit" in record and record["unit"]:
            unit = record["unit"]
        evidence = {}
        cell_type = record.get("Cell type") or record.get("cell_type")
        reliability = record.get("Reliability") or record.get("reliability")
        if cell_type:
            evidence["cell_type"] = cell_type
        if reliability:
            evidence["reliability"] = reliability
        records.append(
            {
                "gene_id": str(gene),
                "tissue": str(tissue),
                "expression": float(level) if level is not None and pd.notna(level) else None,
                "unit": str(unit),
                "source": "HPA",
                "evidence": json.dumps(evidence) if evidence else json.dumps({}),
            }
        )
    result = pd.DataFrame(records, columns=TISSUE_COLUMNS)
    result.drop_duplicates(subset=["gene_id", "tissue", "source"], inplace=True)
    return result


def load_plae(path: Optional[str | Path]) -> pd.DataFrame:
    df = _read_expression(path)
    if df.empty:
        return pd.DataFrame(columns=TISSUE_COLUMNS)

    records: List[Mapping[str, object]] = []
    for record in df.to_dict(orient="records"):
        gene = record.get("gene_id") or record.get("gene")
        tissue = record.get("compartment") or record.get("tissue") or record.get("dataset")
        level = record.get("log2_tpm") or record.get("expression") or record.get("avg_log2")
        if not gene or not tissue:
            continue
        evidence = {}
        dataset = record.get("dataset")
        if dataset:
            evidence["dataset"] = dataset
        cell_type = record.get("cell_type")
        if cell_type:
            evidence["cell_type"] = cell_type
        records.append(
            {
                "gene_id": str(gene),
                "tissue": str(tissue),
                "expression": float(level) if level is not None and pd.notna(level) else None,
                "unit": "log2(TPM)",
                "source": "PLAE",
                "evidence": json.dumps(evidence) if evidence else json.dumps({}),
            }
        )
    result = pd.DataFrame(records, columns=TISSUE_COLUMNS)
    result.drop_duplicates(subset=["gene_id", "tissue", "source"], inplace=True)
    return result


__all__ = ["load_hpa", "load_plae", "TISSUE_COLUMNS"]
