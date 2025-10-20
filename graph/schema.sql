CREATE TABLE IF NOT EXISTS gene (
        gene_id TEXT PRIMARY KEY,
        symbol TEXT NOT NULL,
        ensembl_id TEXT,
        uniprot_id TEXT,
        species TEXT DEFAULT 'Homo sapiens'
);

CREATE TABLE IF NOT EXISTS protein_edge (
        edge_id BIGSERIAL PRIMARY KEY,
        src_gene_id TEXT NOT NULL,
        dst_gene_id TEXT NOT NULL,
        relation TEXT NOT NULL,
        sign TEXT,
        direct BOOLEAN DEFAULT FALSE,
        evidence TEXT,
        source TEXT NOT NULL,
        source_reference TEXT
);

CREATE INDEX IF NOT EXISTS ix_protein_edge_src ON protein_edge(src_gene_id);

CREATE INDEX IF NOT EXISTS ix_protein_edge_dst ON protein_edge(dst_gene_id);

CREATE INDEX IF NOT EXISTS ix_protein_edge_source ON protein_edge(source);

CREATE TABLE IF NOT EXISTS pathway_member (
        pathway_id TEXT NOT NULL,
        gene_id TEXT NOT NULL,
        pathway_name TEXT,
        source TEXT NOT NULL,
        evidence TEXT
);

CREATE INDEX IF NOT EXISTS ix_pathway_member_pathway ON pathway_member(pathway_id);

CREATE INDEX IF NOT EXISTS ix_pathway_member_gene ON pathway_member(gene_id);

CREATE TABLE IF NOT EXISTS disease_gene (
        disease_id TEXT NOT NULL,
        gene_id TEXT NOT NULL,
        evidence_type TEXT NOT NULL,
        score DOUBLE PRECISION,
        source TEXT NOT NULL,
        evidence JSONB
);

CREATE INDEX IF NOT EXISTS ix_disease_gene_disease ON disease_gene(disease_id);

CREATE INDEX IF NOT EXISTS ix_disease_gene_gene ON disease_gene(gene_id);

CREATE TABLE IF NOT EXISTS tissue_expr (
        record_id BIGSERIAL PRIMARY KEY,
        gene_id TEXT NOT NULL,
        tissue TEXT NOT NULL,
        expression DOUBLE PRECISION,
        unit TEXT DEFAULT 'TPM',
        source TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS ix_tissue_expr_gene ON tissue_expr(gene_id);

CREATE INDEX IF NOT EXISTS ix_tissue_expr_tissue ON tissue_expr(tissue);

CREATE TABLE IF NOT EXISTS drug (
        drug_id TEXT PRIMARY KEY,
        preferred_name TEXT NOT NULL,
        synonyms TEXT,
        smiles TEXT,
        inchikey TEXT,
        approval_status TEXT,
        source TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS drug_target (
        drug_id TEXT NOT NULL,
        target_id TEXT NOT NULL,
        target_type TEXT,
        action TEXT,
        affinity DOUBLE PRECISION,
        affinity_unit TEXT,
        moa_category TEXT,
        source TEXT NOT NULL,
        evidence JSONB
);

CREATE INDEX IF NOT EXISTS ix_drug_target_drug ON drug_target(drug_id);

CREATE INDEX IF NOT EXISTS ix_drug_target_target ON drug_target(target_id);

CREATE TABLE IF NOT EXISTS safety_ae (
        drug_id TEXT NOT NULL,
        adverse_event TEXT NOT NULL,
        report_count INTEGER,
        proportional_reporting_ratio DOUBLE PRECISION,
        source TEXT NOT NULL,
        evidence JSONB
);

CREATE TABLE IF NOT EXISTS label_section (
        drug_id TEXT NOT NULL,
        section TEXT NOT NULL,
        text TEXT NOT NULL,
        source TEXT NOT NULL,
        extracted_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE TABLE IF NOT EXISTS trial (
        nct_id TEXT PRIMARY KEY,
        title TEXT NOT NULL,
        status TEXT,
        phase TEXT,
        conditions TEXT,
        enrollment INTEGER,
        interventions TEXT,
        last_updated TIMESTAMPTZ,
        source TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS ix_trial_status ON trial(status);

CREATE INDEX IF NOT EXISTS ix_trial_phase ON trial(phase);
