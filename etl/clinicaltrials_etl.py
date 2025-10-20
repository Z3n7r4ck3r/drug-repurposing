"""Parser for ClinicalTrials.gov v2 exports."""
from __future__ import annotations

from pathlib import Path
from typing import List, Mapping, Optional

import pandas as pd

TRIAL_COLUMNS = ["nct_id", "title", "status", "phase", "conditions", "enrollment", "interventions", "last_updated", "source"]


def load_clinicaltrials(path: Optional[str | Path]) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame(columns=TRIAL_COLUMNS)
    file_path = Path(path)
    if not file_path.exists():
        return pd.DataFrame(columns=TRIAL_COLUMNS)
    if file_path.suffix.lower() in {".json", ".ndjson"}:
        df = pd.read_json(file_path, lines=file_path.suffix.lower() == ".ndjson")
    else:
        sep = "\t" if file_path.suffix.lower() in {".tsv", ".txt"} else ","
        df = pd.read_csv(file_path, sep=sep)
    if df.empty:
        return pd.DataFrame(columns=TRIAL_COLUMNS)

    records: List[Mapping[str, object]] = []
    for record in df.to_dict(orient="records"):
        nct_id = record.get("nct_id") or record.get("NCTId")
        title = record.get("brief_title") or record.get("BriefTitle") or record.get("title")
        if not nct_id or not title:
            continue
        status = record.get("overall_status") or record.get("OverallStatus")
        phase = record.get("phase") or record.get("Phase")
        conditions = record.get("conditions") or record.get("Condition")
        interventions = record.get("interventions") or record.get("InterventionName")
        updated = record.get("last_update_posted_date") or record.get("LastUpdatePostDate")
        enrollment = record.get("enrollment") or record.get("EnrollmentCount")
        try:
            enrollment_value = int(enrollment) if enrollment is not None and str(enrollment).strip() else None
        except ValueError:
            enrollment_value = None
        records.append(
            {
                "nct_id": str(nct_id),
                "title": str(title),
                "status": str(status) if status else "",
                "phase": str(phase) if phase else "",
                "conditions": "; ".join(conditions) if isinstance(conditions, list) else str(conditions) if conditions else "",
                "enrollment": enrollment_value,
                "interventions": "; ".join(interventions) if isinstance(interventions, list) else str(interventions) if interventions else "",
                "last_updated": str(updated) if updated else None,
                "source": "ClinicalTrials.gov",
            }
        )
    result = pd.DataFrame(records, columns=TRIAL_COLUMNS)
    result.drop_duplicates(subset=["nct_id"], inplace=True)
    return result


__all__ = ["load_clinicaltrials", "TRIAL_COLUMNS"]
