"""ClinicalTrials.gov helper utilities."""
from __future__ import annotations

import logging
from typing import Iterable, List

import pandas as pd
import requests

API_URL = "https://clinicaltrials.gov/api/v2/studies"


def fetch_trials(conditions: Iterable[str], recruitment: str = "Recruiting", max_records: int = 200) -> pd.DataFrame:
    """Fetch trial metadata for the supplied conditions."""

    records: List[dict] = []
    for condition in conditions:
        params = {
            "format": "json",
            "cond": condition,
            "fields": "NCTId,BriefTitle,Condition,OverallStatus,StudyType,StartDate,Phase",
            "pageSize": min(max_records, 100),
        }
        try:
            response = requests.get(API_URL, params=params, timeout=60)
            response.raise_for_status()
        except Exception as exc:  # pragma: no cover - network
            logging.warning("Failed to fetch trials for %s: %s", condition, exc)
            continue
        payload = response.json()
        studies = payload.get("studies", [])
        for study in studies:
            attribs = study.get("protocolSection", {}).get("identificationModule", {})
            conditions_list = study.get("protocolSection", {}).get("conditionsModule", {}).get("conditions", [])
            status = study.get("protocolSection", {}).get("statusModule", {}).get("overallStatus")
            if recruitment and status and recruitment.lower() not in status.lower():
                continue
            phases = study.get("protocolSection", {}).get("designModule", {}).get("phases", [])
            if isinstance(phases, list):
                phases = "|".join(phases)
            records.append(
                {
                    "nct_id": attribs.get("nctId"),
                    "title": attribs.get("briefTitle"),
                    "conditions": "|".join(conditions_list),
                    "status": status,
                    "phase": phases,
                }
            )
    return pd.DataFrame(records)


def summarise_trials(trials: pd.DataFrame) -> pd.DataFrame:
    if trials.empty:
        return pd.DataFrame(columns=["condition", "active_trials"])

    exploded = trials.assign(condition_list=trials["conditions"].str.split("|") if "conditions" in trials else [[]])
    exploded = exploded.explode("condition_list")
    exploded["condition_list"].fillna("", inplace=True)
    exploded = exploded[exploded["condition_list"].str.len() > 0]
    return (
        exploded.groupby("condition_list")
        .agg(active_trials=("nct_id", "nunique"))
        .reset_index()
        .rename(columns={"condition_list": "condition"})
    )
