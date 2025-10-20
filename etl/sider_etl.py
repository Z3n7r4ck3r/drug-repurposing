"""Parser for SIDER adverse event data."""
from __future__ import annotations

import json
from pathlib import Path
from typing import List, Mapping, Optional

import pandas as pd

SAFETY_COLUMNS = ["drug_id", "adverse_event", "report_count", "proportional_reporting_ratio", "source", "evidence"]


def load_sider(path: Optional[str | Path]) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame(columns=SAFETY_COLUMNS)
    file_path = Path(path)
    if not file_path.exists():
        return pd.DataFrame(columns=SAFETY_COLUMNS)

    sep = "\t" if file_path.suffix.lower() in {".tsv", ".txt"} else ","
    df = pd.read_csv(file_path, sep=sep)
    if df.empty:
        return pd.DataFrame(columns=SAFETY_COLUMNS)

    records: List[Mapping[str, object]] = []
    for record in df.to_dict(orient="records"):
        drug = record.get("drug_id") or record.get("drug") or record.get("stitch_id_flat")
        event = record.get("adverse_event") or record.get("side_effect_name") or record.get("meddra_concept")
        if not drug or not event:
            continue
        reports = record.get("reports") or record.get("frequency")
        prr = record.get("prr") or record.get("proportional_reporting_ratio")
        evidence = {}
        if "lower_confidence" in record and pd.notna(record["lower_confidence"]):
            evidence["lower_confidence"] = float(record["lower_confidence"])
        if "upper_confidence" in record and pd.notna(record["upper_confidence"]):
            evidence["upper_confidence"] = float(record["upper_confidence"])
        records.append(
            {
                "drug_id": str(drug),
                "adverse_event": str(event),
                "report_count": int(reports) if reports is not None and str(reports).isdigit() else None,
                "proportional_reporting_ratio": float(prr) if prr is not None and pd.notna(prr) else None,
                "source": "SIDER",
                "evidence": json.dumps(evidence) if evidence else json.dumps({}),
            }
        )
    result = pd.DataFrame(records, columns=SAFETY_COLUMNS)
    result.drop_duplicates(subset=["drug_id", "adverse_event"], inplace=True)
    return result


__all__ = ["load_sider", "SAFETY_COLUMNS"]
