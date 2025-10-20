from __future__ import annotations

import json
import re
import shutil
from collections import defaultdict, deque
from datetime import date
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence, Set, Tuple

import pandas as pd

DEFAULT_OUTPUT_ROOT = Path("data/ontology")

DEFAULT_OCULAR_ROOTS = {
    "HP:0000478": "eye",  # Abnormality of the eye
    "HP:0000568": "retina",
    "HP:0000575": "choroid",
    "HP:0000491": "cornea",
    "HP:0000592": "optic_nerve",
    "HP:0012373": "uvea",
    "HP:0000615": "lacrimal_system",
}

def parse_obo(filename, allowed_prefixes):
    """
    Parses an OBO file, yielding only the terms that match the allowed prefixes.
    """
    with open(filename, 'r', encoding='utf-8') as f:
        term = {}
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line == '[Term]':
                if term:
                    # Check if the term's ID has one of the allowed prefixes
                    term_id = term.get('id', [None])[0]
                    if term_id and any(term_id.startswith(p) for p in allowed_prefixes):
                        yield term
                term = {}
            elif ':' in line:
                key, value = line.split(':', 1)
                if key.strip() not in term:
                    term[key.strip()] = []
                term[key.strip()].append(value.strip())
    if term:
        term_id = term.get('id', [None])[0]
        if term_id and any(term_id.startswith(p) for p in allowed_prefixes):
            yield term


def _extract_hpo_relationships(phenotype_terms):
    """Return parent/child maps for the supplied HPO terms."""
    parents = defaultdict(set)
    for term in phenotype_terms:
        term_id = term.get("id", [None])[0]
        if not term_id:
            continue
        for rel in term.get("is_a", []):
            parent_id = rel.split("!")[0].strip()
            if parent_id.startswith("HP"):
                parents[term_id].add(parent_id)

    children = defaultdict(set)
    for child_id, parent_set in parents.items():
        for parent_id in parent_set:
            children[parent_id].add(child_id)

    return parents, children


def _collect_descendants(children: Mapping[str, Set[str]], root_ids: Iterable[str]) -> Set[str]:
    """Breadth-first traversal collecting all descendants for the root IDs."""
    visited = set()
    queue = deque(root_ids)

    while queue:
        current = queue.popleft()
        if current in visited:
            continue
        visited.add(current)
        for child in children.get(current, set()):
            if child not in visited:
                queue.append(child)

    return visited


def build_ocular_descendants(
    children: Mapping[str, Set[str]],
    ocular_roots: Mapping[str, str] | None = None,
) -> Dict[str, Set[str]]:
    """Return the descendant sets for each ocular root."""

    roots = ocular_roots or DEFAULT_OCULAR_ROOTS
    return {root: _collect_descendants(children, {root}) for root in roots}


def annotate_phenotypes(
    phenotype_df: pd.DataFrame,
    ocular_descendants: Mapping[str, Set[str]],
    ocular_roots: Mapping[str, str] | None = None,
) -> Tuple[pd.DataFrame, Sequence[Dict[str, str]]]:
    """Annotate the phenotype dataframe with ocular scope metadata."""

    roots = ocular_roots or DEFAULT_OCULAR_ROOTS

    ocular_phenotypes: Set[str] = set().union(*ocular_descendants.values())

    def infer_scope(hpo_id: str) -> str:
        scopes = [scope for root, scope in roots.items() if hpo_id in ocular_descendants[root]]
        return ",".join(sorted(set(scopes))) if scopes else ""

    phenotype_df = phenotype_df.copy()
    phenotype_df["is_ocular"] = phenotype_df["id"].isin(ocular_phenotypes)
    phenotype_df["ocular_scope"] = phenotype_df["id"].apply(infer_scope)

    subgraph_rows = []
    for root_id, scope in roots.items():
        for phenotype_id in sorted(ocular_descendants.get(root_id, set())):
            subgraph_rows.append(
                {
                    "root_id": root_id,
                    "root_scope": scope,
                    "phenotype_id": phenotype_id,
                }
            )

    return phenotype_df, subgraph_rows


def build_manifest(
    release: str,
    disease_df: pd.DataFrame,
    phenotype_df: pd.DataFrame,
    bridge_df: pd.DataFrame,
    ocular_descendants: Mapping[str, Set[str]],
    ocular_roots: Mapping[str, str] | None = None,
    source_files: Mapping[str, str] | None = None,
) -> Dict[str, object]:
    roots = ocular_roots or DEFAULT_OCULAR_ROOTS
    manifest = {
        "release": release,
        "disease_records": int(len(disease_df)),
        "phenotype_records": int(len(phenotype_df)),
        "ocular_phenotype_records": int(phenotype_df["is_ocular"].sum()),
        "disease_phenotype_records": int(len(bridge_df)),
        "ocular_roots": [],
    }

    if source_files:
        manifest["source_files"] = dict(source_files)

    for root_id, scope in roots.items():
        manifest["ocular_roots"].append(
            {
                "root_id": root_id,
                "scope": scope,
                "descendant_count": int(len(ocular_descendants.get(root_id, set()))),
            }
        )

    return manifest


def process_ontologies(
    ontology_dir: Path | str = Path("."),
    output_root: Path | str = DEFAULT_OUTPUT_ROOT,
    release: str | None = None,
    ensure_latest_copy: bool = True,
):
    """
    Processes the MONDO, DOID, and HPO ontologies to create structured tables
    and a robust disease-phenotype mapping.
    """
    ontology_dir = Path(ontology_dir)
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    release = release or date.today().isoformat()
    output_dir = output_root / release
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Parsing ontology OBO files with strict filtering (release={release})...")

    # --- Create disease table (from MONDO and DOID only) ---
    source_files = {
        "mondo": str((ontology_dir / "mondo.obo").resolve()),
        "doid": str((ontology_dir / "doid.obo").resolve()),
        "hp": str((ontology_dir / "hp.obo").resolve()),
        "hpo_annotations": str((ontology_dir / "phenotype.hpoa").resolve()),
    }

    disease_terms = list(
        parse_obo(ontology_dir / "mondo.obo", ["MONDO", "DOID"])
    ) + list(parse_obo(ontology_dir / "doid.obo", ["MONDO", "DOID"]))

    diseases_data = []
    for term in disease_terms:
        # Additional check to ensure we only process disease terms
        if 'id' in term and (term['id'][0].startswith("MONDO") or term['id'][0].startswith("DOID")):
            synonyms = "|".join([s.split('"')[1] for s in term.get('synonym', []) if '"' in s])

            clean_xrefs = []
            for xref in term.get('xref', []):
                match = re.match(r'([^ ]+)', xref)
                if match:
                    clean_xrefs.append(match.group(1))

            diseases_data.append({
                'id': term['id'][0],
                'name': term['name'][0],
                'synonyms': synonyms,
                'xrefs': "|".join(clean_xrefs)
            })

    disease_df = pd.DataFrame(diseases_data)
    disease_df = disease_df.drop_duplicates(subset='id', keep='first')
    disease_path = output_dir / "disease.csv"
    disease_df.to_csv(disease_path, index=False)
    print(f"Created {disease_path} with {len(disease_df):,} records.")

    # --- Create phenotype table (from HPO only) ---
    phenotype_terms = list(parse_obo(ontology_dir / "hp.obo", ["HP"]))
    phenotypes_data = []
    for term in phenotype_terms:
        if 'id' in term and term['id'][0].startswith("HP"):
            phenotypes_data.append({
                'id': term['id'][0],
                'name': term['name'][0]
            })

    parents, children = _extract_hpo_relationships(phenotype_terms)
    ocular_descendants = build_ocular_descendants(children)

    phenotype_df = pd.DataFrame(phenotypes_data)
    phenotype_df, subgraph_rows = annotate_phenotypes(phenotype_df, ocular_descendants)

    phenotype_path = output_dir / "phenotype.csv"
    phenotype_df.to_csv(phenotype_path, index=False)
    print(
        "Created %s with %s records (ocular phenotypes: %s)."
        % (phenotype_path, f"{len(phenotype_df):,}", f"{phenotype_df['is_ocular'].sum():,}")
    )

    ocular_pheno_path = output_dir / "ocular_phenotypes.csv"
    phenotype_df.loc[phenotype_df["is_ocular"]].to_csv(ocular_pheno_path, index=False)
    print(f"Saved ocular-specific phenotype annotations to {ocular_pheno_path}.")

    ocular_subgraph_path = output_dir / "ocular_subgraphs.csv"
    pd.DataFrame(subgraph_rows).to_csv(ocular_subgraph_path, index=False)
    print(f"Saved ocular scope subgraph relationships to {ocular_subgraph_path}.")

    # --- Build the external to internal lookup table ---
    print("Building external ID to internal disease ID lookup table...")

    def normalize(ext_id):
        if not isinstance(ext_id, str): return None
        ext_id = ext_id.strip().replace("OMIMPS:", "OMIM:")
        if ext_id.startswith("OMIM:"):
            parts = ext_id.split(':')
            if len(parts) > 1 and parts[1].isdigit():
                return f"OMIM:{parts[1].zfill(6)}"
        if ext_id.startswith("ORPHA:"):
            parts = ext_id.split(':')
            if len(parts) > 1 and parts[1].isdigit():
                return f"ORPHA:{int(parts[1])}"
        return ext_id

    ext2disease = {}
    for row in disease_df.itertuples():
        key = normalize(row.id)
        if key not in ext2disease: ext2disease[key] = row.id
        if pd.notna(row.xrefs):
            for ext in row.xrefs.split('|'):
                key = normalize(ext)
                if key and key not in ext2disease: ext2disease[key] = row.id

    print(f"Lookup table created with {len(ext2disease):,} entries.")

    # --- Analyze phenotype.hpoa using the lookup dictionary ---
    print("Processing phenotype.hpoa to create disease-phenotype bridge...")

    hpoa = pd.read_csv(
        ontology_dir / "phenotype.hpoa", sep="\t", skiprows=4, header=0,
        names=["database_id", "disease_name", "qualifier", "hpo_id", "reference", "evidence", "onset", "frequency", "sex", "modifier", "aspect", "biocuration"],
        dtype=str
    )[lambda x: (x["qualifier"] != "NOT") & x["hpo_id"].notna()].copy()

    hpoa["disease_id"] = hpoa["database_id"].apply(lambda db_id: ext2disease.get(normalize(db_id)))
    hpoa["phenotype_id"] = hpoa["hpo_id"].str.strip()

    bridge = hpoa.dropna(subset=["disease_id"])[["disease_id", "phenotype_id"]].drop_duplicates()
    bridge_path = output_dir / "disease_phenotype.csv"
    bridge.to_csv(bridge_path, index=False)
    print(f"Created {bridge_path} with {len(bridge):,} records.")

    manifest = build_manifest(
        release,
        disease_df,
        phenotype_df,
        bridge,
        ocular_descendants,
        source_files=source_files,
    )

    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote ontology manifest to {manifest_path}.")

    if ensure_latest_copy:
        latest_dir = output_root / "latest"
        if latest_dir.exists():
            shutil.rmtree(latest_dir)
        shutil.copytree(output_dir, latest_dir)
        print(f"Updated latest ontology snapshot at {latest_dir}.")

if __name__ == "__main__":
    process_ontologies()
