"""Streamlit front-end to explore repurposing outputs."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

SEED_PATH = Path("disease_seeds.csv")
SCORE_PATH = Path("target_scores.csv")
DRUG_TARGET_PATH = Path("drug_target_mapping.csv")


@st.cache_data
def load_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def main() -> None:
    st.title("Rx Repurpose Dashboard")

    seeds = load_table(SEED_PATH)
    scores = load_table(SCORE_PATH)
    targets = load_table(DRUG_TARGET_PATH)

    disease_list = sorted(seeds["disease_id"].unique()) if not seeds.empty else []
    selected_disease = st.selectbox("Select disease", disease_list)

    if selected_disease:
        disease_seeds = seeds[seeds["disease_id"] == selected_disease]
        st.subheader("Seed genes")
        st.dataframe(disease_seeds[["gene_symbol", "score", "source"]])

        if not scores.empty:
            if "disease_id" in scores.columns:
                disease_scores = scores[scores["disease_id"].astype(str) == selected_disease]
            elif "seed_set" in scores.columns:
                disease_scores = scores[scores["seed_set"] == selected_disease]
            else:
                disease_scores = scores
            st.subheader("Propagated targets")
            st.dataframe(disease_scores.head(50))

        if not targets.empty:
            linked = targets[targets["disease_id"] == selected_disease] if "disease_id" in targets.columns else targets
            st.subheader("Drug-targets")
            st.dataframe(linked.head(50))


if __name__ == "__main__":
    main()
