# Remaining Steps to Achieve a SOTA Ocular Drug Repurposing Platform

This roadmap captures the concrete work still required to fulfil `plan.txt` and deliver a state-of-the-art pipeline for ocular drug repurposing. Steps are ordered to respect upstream dependencies; each section lists the key datasets, implementation artefacts, validation gates, and deliverables needed before advancing.

## 1. Ontology Scope Hardening
- **Data**: Refresh MONDO, DOID, HPO dumps; incorporate SNOMED CT ocular subsets if license permits.
- **Implementation**: Extend `etl/ontology_parser.py` to materialise explicit ocular subgraphs (retina, cornea, uvea, optic nerve, lacrimal) with persistent identifiers; unify synonyms and cross-references.
- **Validation**: Automated acceptance queries for AMD, glaucoma, uveitis, retinal vein occlusion, keratitis returning MONDO/DOID plus ocular HPO phenotypes.
- **Deliverables**: Versioned ontology parquet/CSV exports; regression tests locking expected ocular coverage.

## 2. Data Lake & Knowledge Graph Completion
- **Database**: Stand up PostgreSQL schema (`gene`, `protein_edge`, `pathway_member`, `disease_gene`, `tissue_expr`, `drug`, `drug_target`, `safety_ae`, `label_section`, `trial`) and migration scripts.
- **ETL Coverage**:
  - Causal networks: OmniPath, SIGNOR (activation/inhibition, direct flag).
  - PPI: STRING (physical score tiers).
  - Pathways: Reactome gene sets.
  - Genetics: OpenTargets, GWAS Catalog, DisGeNET.
  - Expression: HPA Eye/Retina, PLAE/scEiaD bulk & single-cell.
  - Pharmacology: ChEMBL, DrugCentral, IUPHAR, RxNorm/RxClass harmonisation.
  - Safety/Labels: FAERS, DailyMed SPL, SIDER.
  - Trials: ClinicalTrials.gov v2.
- **Implementation**: Modular loaders under `etl/` with ID normalization (HGNC/Ensembl/UniProt, RxNorm/ChEMBL) and staging-to-core pipelines; orchestrate via Airflow/Snakemake DAG.
- **Validation**: `protein_edge` signed direct edges >200k; cross-source drug_target coverage benchmarks; automated schema/data quality tests.
- **Deliverables**: Documented ETL DAG, data versioning metadata (release dates, checksums), Dockerised Postgres snapshot.

## 3. Disease Module & Expression Signatures
- **Status**: Core differential expression + DerSimonian–Laird meta-analysis pipeline and CLI delivered; remaining work centres on full data ingestion, retina-specialised weighting, and persistence into the shared store.
- **Datasets**: GEO/ArrayExpress ocular studies, retina-specific PLAE datasets, OpenTargets and GWAS summary statistics (ingestion + cataloguing outstanding).
- **Implementation**: Automate study discovery/selection, wire GEO/PLAE loaders into the new pipeline, add retina/cell-type weighting, and write normalised outputs into Postgres `disease_gene` tables alongside genetics and curated seeds.
- **Validation**: Produce QC reports (batch effects, cell-type coverage), benchmark against known ocular gene sets, and extend unit/integration tests to cover end-to-end data ingestion plus score persistence.
- **Deliverables**: Reproducible notebooks or Snakemake rules for study sourcing, cached signature artefacts with manifests, and documentation detailing dataset selection and QC outcomes.

## 4. Network Scoring Enhancements
- **Algorithms**: Finalise propagation, random walk with restart, interconnectivity modules; add alternative propagation baselines (heat kernel, PageRank) for comparison.
- **Implementation**: Optimise for large-scale graphs (sparse matrices, GPU optional); integrate with disease modules to output ranked targets with provenance.
- **Validation**: Hold-out evaluation using known approved therapies, sensitivity analyses vs. parameter sweeps, unit tests for numerical stability.
- **Deliverables**: `scoring/` CLI producing `target_score` tables with metadata, reproducible benchmarking reports.

## 5. Drug Target Mapping & Directionality
- **Data Harmonisation**: Merge ChEMBL, DrugCentral, IUPHAR, RxNorm; resolve mechanism-of-action, action type, affinity metrics.
- **Implementation**: Extend `mapping/drug_target_map.py` to annotate direction of effect (agonist/antagonist/etc.), molecular modality, and clinical phase.
- **Validation**: Consistency checks against gold-standard lists (e.g., DrugBank sample), coverage statistics per mechanism, automated tests for ID reconciliation.
- **Deliverables**: Unified `drug_target` table, documented mapping rules, change logs for data refreshes.

## 6. Signature Reversal & Perturbagen Integration
- **Datasets**: LINCS L1000 (Level 5), iLINCS contrasts, DSigDB gene sets.
- **Implementation**: Normalise perturbational signatures, compute connectivity scores (tau/cmap), aggregate by drug and cell type; expose CLI to produce `signature_score` outputs.
- **Validation**: Replicate known positive controls (e.g., approved AMD therapies), cross-validate against disease signatures.
- **Deliverables**: Processed signature repository, provenance manifest, tests covering score calculation.

## 7. Developability & Safety Analytics
- **Pipelines**: Robust FAERS signal detection, DailyMed NER with guard-rails (RxNorm/HPO constraints), SIDER adverse event parsing, ClinicalTrials.gov extraction with ocular filtering.
- **Implementation**: Compose `developability/` modules into scoring functions producing `D_d`; integrate trial status, dosage forms, route of administration.
- **Validation**: Manual spot-check dashboards, statistical thresholds for adverse event enrichment, regression tests on parsers.
- **Deliverables**: Safety score datasets, documentation of guard-rails, alerting for data ingestion drift.

## 8. Evidence Integration & Modeling
- **Formula**: Implement final ranking `F_d = w1·R_d + w2·D_d + w3·tau` with learnable or configurable weights; support ensemble blending and calibration.
- **Implementation**: `model/` modules for feature engineering, training (e.g., gradient boosting, logistic regression), cross-validation on historical successes.
- **Validation**: Backtesting on ocular indications with known repurposed drugs, ablation studies, uncertainty quantification.
- **Deliverables**: Reproducible training scripts, model cards, evaluation reports.

## 9. API & UI Completion
- **Backend**: Expand FastAPI to cover search, graph traversal, evidence cards, provenance; add authentication and rate limiting if exposed externally.
- **Frontend**: Rich Streamlit (or alternative) dashboards with network visualisations, drill-down evidence, signature comparisons.
- **Validation**: End-to-end integration tests, usability walkthroughs, screenshot baselines, accessibility checks.
- **Deliverables**: Deployment-ready Docker images, CI pipeline for automated testing, documentation for operators.

## 10. Governance, QA, and DevOps
- **Versioning**: Introduce DVC or similar for dataset snapshots; maintain release manifest with source versions, hashes, and licensing notes.
- **Testing**: Extend unit/integration coverage, add data integrity tests, bootstrapped ranking stability checks.
- **Ops**: Harden Airflow/Snakemake orchestration, add monitoring/logging, implement secure credential handling.
- **Deliverables**: QA checklist, incident response playbook, CI/CD workflows.

## 11. Documentation & Publication Readiness
- **Docs**: Comprehensive README, technical whitepaper, API docs (OpenAPI/Swagger), user guides for analysts.
- **Reproducibility**: End-to-end runbook, sample datasets, scripted environment setup (Docker/Conda), licensing compliance summary.
- **Outreach**: Benchmark comparisons vs. existing repurposing tools, highlight SOTA performance metrics, prepare for potential open publication.

---

**Execution Guidance**
1. Treat each section as a milestone with its own feature branch, tests, and review checklist.
2. Maintain a changelog mapping commits to roadmap items and data releases.
3. Prioritise data quality and provenance; every score must be traceable to source evidence.
4. Continuously benchmark against industry/academic baselines to validate SOTA claims.
