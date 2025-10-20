"""Structured extraction of ocular safety snippets from DailyMed SPL."""
from __future__ import annotations

import re
from typing import Dict, Iterable, List

import pandas as pd

OCULAR_PATTERNS = [
    re.compile(pattern, re.IGNORECASE)
    for pattern in [
        r"ocular",
        r"ophthalm",
        r"eye",
        r"retin",
        r"cornea",
        r"uveitis",
        r"glaucoma",
    ]
]


def load_dailymed_sections(path: str = "dailymed_processed_data.csv") -> pd.DataFrame:
    try:
        return pd.read_csv(path)
    except FileNotFoundError:
        return pd.DataFrame()


def extract_ocular_snippets(sections: pd.DataFrame, patterns: Iterable[re.Pattern] = OCULAR_PATTERNS) -> pd.DataFrame:
    if sections.empty:
        return pd.DataFrame(columns=["drug_name", "section", "snippet"])

    records: List[Dict[str, str]] = []
    for row in sections.to_dict(orient="records"):
        for section, text in row.items():
            if section.lower() in {"adverse reactions", "warnings and precautions", "dosage and administration"}:
                if not isinstance(text, str):
                    continue
                for pattern in patterns:
                    if pattern.search(text):
                        snippet = pattern.search(text)
                        start = max(snippet.start() - 80, 0)
                        end = snippet.end() + 120
                        records.append(
                            {
                                "drug_name": row.get("drug_name", row.get("setid", "unknown")),
                                "section": section,
                                "snippet": text[start:end].strip(),
                            }
                        )
                        break
    return pd.DataFrame(records)
