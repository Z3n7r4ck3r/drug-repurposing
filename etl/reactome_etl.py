"""Utilities to parse Reactome pathway membership exports."""
from __future__ import annotations

from pathlib import Path
from typing import List, Mapping, Optional

import pandas as pd

PATHWAY_COLUMNS = ["pathway_id", "pathway_name", "gene_id", "evidence"]


def _coerce_path(path: Optional[str | Path]) -> Optional[Path]:
    if path is None:
        return None
    candidate = Path(path)
    if candidate.exists():
        return candidate
    return None


def load_reactome(path: Optional[str | Path]) -> pd.DataFrame:
    """Load a Reactome pathway membership file.

    Parameters
    ----------
    path:
        Location of the TSV/CSV file. The export should contain columns similar
        to the Reactome `Ensembl2Reactome_All_Levels` files. If the file is
        absent the function returns an empty data frame.
    """

    file_path = _coerce_path(path)
    if file_path is None:
        return pd.DataFrame(columns=PATHWAY_COLUMNS)

    sep = "\t" if file_path.suffix.lower() in {".tsv", ".txt"} else ","
    df = pd.read_csv(file_path, sep=sep)
    if df.empty:
        return pd.DataFrame(columns=PATHWAY_COLUMNS)

    records: List[Mapping[str, str]] = []
    for record in df.to_dict(orient="records"):
        pathway_id = record.get("pathway_id") or record.get("Pathway identifier") or record.get("Pathway stId")
        if not pathway_id:
            for key in ("Pathway Identifier", "stId"):
                value = record.get(key)
                if value:
                    pathway_id = value
                    break
        gene = record.get("gene_id") or record.get("Entity identifier") or record.get("Entity Identifier")
        if not gene:
            for key in ("Ensembl identifier", "Ensembl Identifier", "Entity", "Entity ID"):
                value = record.get(key)
                if value:
                    gene = value
                    break
        if not pathway_id or not gene:
            continue

        pathway_name = (
            record.get("pathway_name")
            or record.get("Pathway name")
            or record.get("Pathway Name")
            or record.get("Pathway")
        )
        evidence = record.get("Evidence") or record.get("evidence")

        records.append(
            {
                "pathway_id": str(pathway_id),
                "pathway_name": str(pathway_name) if pathway_name else "",
                "gene_id": str(gene),
                "evidence": str(evidence) if evidence else "",
            }
        )

    if not records:
        return pd.DataFrame(columns=PATHWAY_COLUMNS)

    result = pd.DataFrame(records, columns=PATHWAY_COLUMNS)
    result.drop_duplicates(subset=["pathway_id", "gene_id"], inplace=True)
    result["pathway_name"].fillna("", inplace=True)
    result["evidence"].fillna("", inplace=True)
    return result


__all__ = ["load_reactome", "PATHWAY_COLUMNS"]
