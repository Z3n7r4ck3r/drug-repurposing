import argparse
from pathlib import Path
from collections import defaultdict

import pandas as pd

OCULAR_KEYWORDS = [
    "ocular",
    "eye",
    "retina",
    "retinal",
    "macula",
    "macular",
    "cornea",
    "uveitis",
    "uvea",
    "glaucoma",
    "optic nerve",
    "choroid",
    "lacrimal",
]


def _load_tables(input_dir: Path):
    disease_df = pd.read_csv(input_dir / "disease.csv")
    phenotype_df = pd.read_csv(input_dir / "phenotype.csv")
    bridge_df = pd.read_csv(input_dir / "disease_phenotype.csv")
    return disease_df, phenotype_df, bridge_df


def _aggregate_phenotypes(phenotype_df, bridge_df):
    ocular_ids = set(
        phenotype_df.loc[phenotype_df["is_ocular"], "id"]
    )
    label_map = phenotype_df.set_index("id")["name"].to_dict()
    scope_map = phenotype_df.set_index("id")["ocular_scope"].to_dict()

    ocular_bridge = bridge_df[bridge_df["phenotype_id"].isin(ocular_ids)]
    aggregated = defaultdict(lambda: {"ids": [], "labels": [], "scopes": set()})

    for row in ocular_bridge.itertuples(index=False):
        phenotype_id = row.phenotype_id
        disease_id = row.disease_id
        aggregated[disease_id]["ids"].append(phenotype_id)
        if phenotype_id in label_map:
            aggregated[disease_id]["labels"].append(label_map[phenotype_id])
        scope_value = scope_map.get(phenotype_id, "")
        if scope_value:
            for scope in scope_value.split(","):
                if scope:
                    aggregated[disease_id]["scopes"].add(scope)

    return aggregated


def tag_ocular_diseases(input_dir: Path, output_dir: Path | None = None):
    """Tag diseases that belong to the ocular scope using ontology relations."""

    output_dir = output_dir or input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    disease_df, phenotype_df, bridge_df = _load_tables(input_dir)
    aggregated = _aggregate_phenotypes(phenotype_df, bridge_df)

    disease_df["is_ocular"] = False
    disease_df["ocular_phenotype_ids"] = ""
    disease_df["ocular_phenotype_labels"] = ""
    disease_df["ocular_scopes"] = ""
    disease_df["tag_source"] = ""

    for idx, row in disease_df.iterrows():
        disease_id = row["id"]
        sources = []
        if disease_id in aggregated:
            payload = aggregated[disease_id]
            unique_ids = sorted(set(payload["ids"]))
            unique_labels = sorted(set(payload["labels"]))
            scopes = sorted(payload["scopes"])

            disease_df.at[idx, "is_ocular"] = True
            disease_df.at[idx, "ocular_phenotype_ids"] = "|".join(unique_ids)
            disease_df.at[idx, "ocular_phenotype_labels"] = "|".join(unique_labels)
            disease_df.at[idx, "ocular_scopes"] = "|".join(scopes)
            sources.append("phenotype")

        if not disease_df.at[idx, "is_ocular"]:
            name = str(row["name"]).lower()
            synonyms = str(row.get("synonyms", "")).lower()
            for keyword in OCULAR_KEYWORDS:
                if keyword in name or keyword in synonyms:
                    disease_df.at[idx, "is_ocular"] = True
                    sources.append("keyword")
                    break
        elif str(row.get("synonyms", "")).strip():
            synonyms = str(row.get("synonyms", "")).lower()
            name = str(row["name"]).lower()
            if any(keyword in name or keyword in synonyms for keyword in OCULAR_KEYWORDS):
                sources.append("keyword")

        if sources:
            disease_df.at[idx, "tag_source"] = "|".join(sorted(set(sources)))

    tagged_path = output_dir / "disease_tagged.csv"
    disease_df.to_csv(tagged_path, index=False)
    print(f"Tagged {int(disease_df['is_ocular'].sum()):,} ocular diseases.")
    print(f"Saved the tagged data to {tagged_path}")

    ocular_summary = disease_df[disease_df["is_ocular"]].copy()
    ocular_summary[[
        "id",
        "name",
        "ocular_scopes",
        "ocular_phenotype_ids",
        "ocular_phenotype_labels",
        "tag_source",
    ]].to_csv(output_dir / "ocular_disease_summary.csv", index=False)
    print(f"Saved ocular disease summary to {output_dir / 'ocular_disease_summary.csv'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("data/ontology/latest"),
        help="Directory containing ontology exports (disease.csv, phenotype.csv, disease_phenotype.csv).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Destination directory for tagged outputs. Defaults to --input-dir.",
    )
    args = parser.parse_args()

    input_dir = args.input_dir
    if not (input_dir / "disease.csv").exists():
        raise FileNotFoundError(
            f"Could not find ontology exports in {input_dir}. Run ontology_parser first."
        )

    tag_ocular_diseases(input_dir, args.output_dir)


if __name__ == "__main__":
    main()
