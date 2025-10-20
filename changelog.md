# Changelog

This log captures the incremental progress toward the SOTA ocular drug repurposing platform defined in `plan.txt` and detailed in `remaining_steps.md`.

## [Unreleased]
### Added
- Versioned ontology export pipeline with ocular scope annotations, subgraph tables, and provenance manifests.
- Snapshot-aware ocular disease tagging and query CLIs aligned with the versioned ontology outputs.
- Declarative PostgreSQL knowledge-graph schema generator with optional apply step.
- Expanded knowledge-graph builder covering causal, pathway, genetics, expression, pharmacology, safety, and trial data feeds with SQLite/PostgreSQL loaders.
- Dataset-specific ETL loaders for Reactome, genetics, expression, RxNorm, SIDER, and ClinicalTrials inputs.
- Airflow orchestration scaffolding and regression tests for schema generation and graph assembly flows.
- Disease-module engine with differential expression meta-analysis, evidence integration, and scoring CLI.

### Progress
- Disease module milestone: core meta-analysis pipeline and CLI are operational; next actions focus on ingesting retina-specific GEO/PLAE cohorts, expanding QC reporting, and wiring outputs into the central `disease_gene` store.

### Planned
Refer to `remaining_steps.md` for the authoritative list of outstanding milestones and acceptance criteria across ontology hardening, data lake completion, disease modules, scoring, signature reversal, developability analytics, evidence integration, API/UI, and governance workstreams.

