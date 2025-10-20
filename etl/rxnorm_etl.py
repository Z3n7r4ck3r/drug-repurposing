"""Helpers for RxNorm / RxClass harmonisation."""
from __future__ import annotations

from pathlib import Path
from typing import List, Mapping, Optional

import pandas as pd

DRUG_COLUMNS = ["drug_id", "preferred_name", "synonyms", "source"]


def _read_table(path: Optional[str | Path]) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame()
    file_path = Path(path)
    if not file_path.exists():
        return pd.DataFrame()
    sep = "\t" if file_path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(file_path, sep=sep)


def load_rxnorm_drugs(path: Optional[str | Path]) -> pd.DataFrame:
    df = _read_table(path)
    if df.empty:
        return pd.DataFrame(columns=DRUG_COLUMNS)

    records: List[Mapping[str, object]] = []
    for record in df.to_dict(orient="records"):
        rxcui = record.get("RXCUI") or record.get("rxcui") or record.get("drug_id")
        name = record.get("STR") or record.get("name") or record.get("preferred_name")
        if not rxcui or not name:
            continue
        synonyms = record.get("synonym") or record.get("SYNONYMS") or record.get("synonyms")
        if isinstance(synonyms, list):
            synonyms_text = ";".join(str(item) for item in synonyms)
        else:
            synonyms_text = str(synonyms) if synonyms else ""
        records.append(
            {
                "drug_id": str(rxcui),
                "preferred_name": str(name),
                "synonyms": synonyms_text,
                "source": "RxNorm",
            }
        )
    result = pd.DataFrame(records, columns=DRUG_COLUMNS)
    result.drop_duplicates(subset=["drug_id"], inplace=True)
    return result


__all__ = ["load_rxnorm_drugs", "DRUG_COLUMNS"]
