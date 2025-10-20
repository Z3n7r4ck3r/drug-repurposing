"""Analyse FAERS/OpenFDA adverse event data for ocular safety signals."""
from __future__ import annotations

import logging
from typing import Iterable

import pandas as pd

OCULAR_TERMS = {
    "retinal detachment",
    "visual acuity reduced",
    "eye pain",
    "ocular hyperaemia",
    "photophobia",
    "blindness",
    "vision blurred",
    "conjunctival haemorrhage",
}


def load_openfda_events(term: str) -> pd.DataFrame:
    try:
        from etl.openfda_etl import get_openfda_adverse_events  # type: ignore
    except Exception as exc:  # pragma: no cover - optional dependency
        logging.warning("openFDA connector unavailable: %s", exc)
        return pd.DataFrame(columns=["safetyreportid", "reactions", "drugs"])
    return get_openfda_adverse_events(term)


def extract_ocular_events(events: pd.DataFrame, ocular_terms: Iterable[str] = OCULAR_TERMS) -> pd.DataFrame:
    if events.empty:
        return pd.DataFrame(columns=["safetyreportid", "reaction", "drug"])

    rows = []
    ocular_terms_lc = {term.lower() for term in ocular_terms}
    for record in events.to_dict(orient="records"):
        reactions = record.get("reactions", []) or []
        drugs = record.get("drugs", []) or []
        for reaction in reactions:
            if reaction and reaction.lower() in ocular_terms_lc:
                for drug in drugs:
                    rows.append({
                        "safetyreportid": record.get("safetyreportid"),
                        "reaction": reaction,
                        "drug": drug,
                    })
    return pd.DataFrame(rows)


def summarise_ocular_events(ocular_events: pd.DataFrame) -> pd.DataFrame:
    if ocular_events.empty:
        return pd.DataFrame(columns=["drug", "reports", "unique_reactions"])
    return (
        ocular_events.groupby("drug")
        .agg(reports=("safetyreportid", "nunique"), unique_reactions=("reaction", lambda x: "|".join(sorted(set(x)))))
        .reset_index()
    )
