"""Utility helpers to query the ocular disease scope derived from ontologies."""
import argparse
from pathlib import Path
from typing import Dict, Iterable, List

import pandas as pd

DEFAULT_QUERIES = [
    "AMD",
    "glaucoma",
    "uveite",
    "occlusione venosa retinica",
    "cheratite",
]


def _require_files(paths: Iterable[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "The following prerequisite files are missing: " + ", ".join(missing)
        )


def _load_tables(base_dir: Path) -> Dict[str, pd.DataFrame]:
    files = {
        "disease": base_dir / "disease.csv",
        "phenotype": base_dir / "phenotype.csv",
        "bridge": base_dir / "disease_phenotype.csv",
        "tagged": base_dir / "disease_tagged.csv",
    }
    _require_files(files.values())

    return {name: pd.read_csv(path) for name, path in files.items()}


def _build_hpo_maps(phenotype_df: pd.DataFrame, bridge_df: pd.DataFrame):
    label_map = phenotype_df.set_index("id")["name"].to_dict()
    ocular_flag = phenotype_df.set_index("id")["is_ocular"].to_dict()
    scope_map = phenotype_df.set_index("id")["ocular_scope"].to_dict()

    disease_to_hpo: Dict[str, List[str]] = {}
    for row in bridge_df.itertuples(index=False):
        disease_to_hpo.setdefault(row.disease_id, []).append(row.phenotype_id)

    return disease_to_hpo, label_map, ocular_flag, scope_map


def _format_hpo_list(hpo_ids: Iterable[str], label_map, ocular_flag, scope_map):
    ocular_entries, other_entries = [], []
    for hpo_id in hpo_ids:
        label = label_map.get(hpo_id, "")
        scope = scope_map.get(hpo_id, "")
        entry = f"{hpo_id} :: {label}" if label else hpo_id
        if scope:
            entry += f" [scope={scope}]"
        if ocular_flag.get(hpo_id, False):
            ocular_entries.append(entry)
        else:
            other_entries.append(entry)
    return ocular_entries, other_entries


def query_diseases(term: str, tables: Dict[str, pd.DataFrame], limit: int = 5) -> pd.DataFrame:
    disease_df = tables["disease"].copy()
    if "synonyms" not in disease_df:
        disease_df["synonyms"] = ""
    disease_df["synonyms"].fillna("", inplace=True)

    mask = (
        disease_df["name"].str.contains(term, case=False, na=False)
        | disease_df["synonyms"].str.contains(term, case=False, na=False)
    )
    result = disease_df[mask].copy()

    tagged_cols = tables["tagged"][
        ["id", "is_ocular", "ocular_scopes", "tag_source"]
    ].copy()
    result = result.merge(tagged_cols, on="id", how="left")

    return result.head(limit)


def run_queries(queries: Iterable[str], limit: int = 5, base_dir: Path | None = None) -> None:
    base_dir = base_dir or Path("data/ontology/latest")
    tables = _load_tables(base_dir)
    disease_to_hpo, label_map, ocular_flag, scope_map = _build_hpo_maps(
        tables["phenotype"], tables["bridge"]
    )

    tagged_info = tables["tagged"].set_index("id")

    for query in queries:
        print(f"\n=== Query: {query} ===")
        matches = query_diseases(query, tables, limit)
        if matches.empty:
            print("No matches found.")
            continue

        for row in matches.itertuples():
            disease_id = row.id
            hpo_ids = disease_to_hpo.get(disease_id, [])
            ocular_entries, other_entries = _format_hpo_list(
                hpo_ids, label_map, ocular_flag, scope_map
            )

            tag_row = tagged_info.loc[disease_id] if disease_id in tagged_info.index else None
            is_ocular = bool(tag_row["is_ocular"]) if tag_row is not None else False
            scopes = tag_row.get("ocular_scopes", "") if tag_row is not None else ""
            sources = tag_row.get("tag_source", "") if tag_row is not None else ""

            print(f"{disease_id} :: {row.name}")
            print(f"  Ocular scope: {is_ocular} (scopes={scopes or '-'}, sources={sources or '-'})")
            if ocular_entries:
                print("  Ocular HPO phenotypes:")
                for item in ocular_entries:
                    print(f"    - {item}")
            if other_entries:
                print("  Additional HPO phenotypes:")
                for item in other_entries:
                    print(f"    - {item}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--disease",
        "-d",
        action="append",
        dest="diseases",
        help="Disease query (can be repeated). Default queries cover the acceptance checks.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=5,
        help="Maximum number of matches per query.",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("data/ontology/latest"),
        help="Directory containing ontology exports and disease tagging outputs.",
    )
    args = parser.parse_args()

    queries = args.diseases if args.diseases else DEFAULT_QUERIES
    run_queries(queries, limit=args.limit, base_dir=args.input_dir)


if __name__ == "__main__":
    main()
