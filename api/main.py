from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd
from fastapi import FastAPI, HTTPException

app = FastAPI(title="Rx Repurpose API", version="0.1.0")

SEED_PATH = Path("disease_seeds.csv")
SCORE_PATH = Path("target_scores.csv")
DRUG_TARGET_PATH = Path("drug_target_mapping.csv")


def _load_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


@app.get("/health")
def health() -> dict:
    return {"status": "ok"}


@app.get("/diseases/{disease_id}")
def disease_detail(disease_id: str) -> dict:
    seeds = _load_csv(SEED_PATH)
    if seeds.empty:
        raise HTTPException(status_code=404, detail="Seed table not found")

    subset = seeds[seeds["disease_id"].str.lower() == disease_id.lower()]
    if subset.empty:
        raise HTTPException(status_code=404, detail="Disease not found")

    return {
        "disease_id": disease_id,
        "seed_genes": subset[["gene_symbol", "score", "source"]].to_dict(orient="records"),
    }


@app.get("/targets/{target_id}")
def target_detail(target_id: str) -> dict:
    scores = _load_csv(SCORE_PATH)
    targets = _load_csv(DRUG_TARGET_PATH)

    score_row: Optional[pd.Series] = None
    if not scores.empty and "target_id" in scores.columns:
        matched = scores[scores["target_id"].str.lower() == target_id.lower()]
        if not matched.empty:
            score_row = matched.iloc[0]

    linked_drugs = []
    if not targets.empty and "target_id" in targets.columns:
        linked_drugs = targets[targets["target_id"].str.lower() == target_id.lower()].to_dict(orient="records")

    if score_row is None and not linked_drugs:
        raise HTTPException(status_code=404, detail="Target not found")

    payload = {"target_id": target_id, "drugs": linked_drugs}
    if score_row is not None and "score" in score_row:
        payload["score"] = float(score_row["score"])
    return payload
