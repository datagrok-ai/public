-- Migration: Add sequence regions, targets, sequence properties/liabilities, and expanded assay types

-- 1. Sequence regions (annotated domains of antibody sequences)
CREATE TABLE biologics.sequence_regions (
  id SERIAL PRIMARY KEY,
  sequence_id INTEGER NOT NULL REFERENCES biologics.sequences(id) ON UPDATE CASCADE ON DELETE CASCADE,
  region_type TEXT NOT NULL,
  start_pos INTEGER NOT NULL CHECK (start_pos >= 1),
  end_pos INTEGER NOT NULL CHECK (end_pos >= start_pos),
  notes TEXT
);

CREATE INDEX ON biologics.sequence_regions (sequence_id);
CREATE INDEX ON biologics.sequence_regions (region_type);

-- 2. Targets (molecular/biological targets for assays)
CREATE TABLE biologics.targets (
  id SERIAL PRIMARY KEY,
  identifier TEXT GENERATED ALWAYS AS ('GROKTGT-' || lpad(id::text, 6, '0')) STORED UNIQUE,
  name TEXT NOT NULL,
  description TEXT
);

-- 3. Add target_id to assay_results
ALTER TABLE biologics.assay_results
  ADD COLUMN IF NOT EXISTS target_id INTEGER REFERENCES biologics.targets(id) ON UPDATE CASCADE ON DELETE SET NULL;

CREATE INDEX IF NOT EXISTS biologics_assay_results_target_id_idx ON biologics.assay_results (target_id);

-- 4. Add category and fit_model to assay_types
ALTER TABLE biologics.assay_types
  ADD COLUMN IF NOT EXISTS category TEXT,
  ADD COLUMN IF NOT EXISTS fit_model TEXT;

-- 5. Add sequence_id to assay_results for direct antibody-level assays
ALTER TABLE biologics.assay_results
  ADD COLUMN IF NOT EXISTS sequence_id INTEGER REFERENCES biologics.sequences(id) ON UPDATE CASCADE ON DELETE SET NULL;

CREATE INDEX IF NOT EXISTS biologics_assay_results_sequence_id_idx ON biologics.assay_results (sequence_id);

-- 6. Sequence properties (computed biophysical properties for hit optimization)
CREATE TABLE biologics.sequence_properties (
  id SERIAL PRIMARY KEY,
  sequence_id INTEGER NOT NULL REFERENCES biologics.sequences(id) ON UPDATE CASCADE ON DELETE CASCADE,
  molecular_weight REAL,
  isoelectric_point REAL,
  hydrophobicity_index REAL,
  charge_at_ph7 REAL,
  aggregation_propensity REAL,
  computed_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX ON biologics.sequence_properties (sequence_id);

-- 7. Sequence liabilities (PTM and developability risk sites)
CREATE TABLE biologics.sequence_liabilities (
  id SERIAL PRIMARY KEY,
  sequence_id INTEGER NOT NULL REFERENCES biologics.sequences(id) ON UPDATE CASCADE ON DELETE CASCADE,
  liability_type TEXT NOT NULL,
  position INTEGER NOT NULL CHECK (position >= 1),
  motif TEXT,
  risk_level TEXT CHECK (risk_level IN ('Low', 'Medium', 'High'))
);

CREATE INDEX ON biologics.sequence_liabilities (sequence_id);
CREATE INDEX ON biologics.sequence_liabilities (liability_type);

-- Sample data: molecular targets
INSERT INTO biologics.targets (name, description) VALUES
  ('HER2', 'Human epidermal growth factor receptor 2'),
  ('PD-1', 'Programmed cell death protein 1'),
  ('PD-L1', 'Programmed death-ligand 1'),
  ('EGFR', 'Epidermal growth factor receptor'),
  ('TNF-alpha', 'Tumor necrosis factor alpha'),
  ('VEGF', 'Vascular endothelial growth factor'),
  ('CD20', 'B-lymphocyte antigen CD20'),
  ('IL-6', 'Interleukin 6'),
  ('CD19', 'B-lymphocyte antigen CD19'),
  ('CTLA-4', 'Cytotoxic T-lymphocyte-associated protein 4');

-- Update existing assay types with categories
UPDATE biologics.assay_types SET category = 'Legacy', fit_model = NULL WHERE category IS NULL;

-- New assay types: comprehensive antibody characterization panel
INSERT INTO biologics.assay_types (name, category, fit_model, description) VALUES
  -- Primary Screening
  ('ELISA Binding OD450', 'Primary Screening', '4PL logistic', 'ELISA optical density at 450 nm'),
  ('ELISA Binding EC50', 'Primary Screening', '4PL logistic', 'ELISA half-maximal effective concentration'),
  ('Flow Cytometry % Positive', 'Primary Screening', '4PL logistic', 'Percentage of positive cells by flow cytometry'),
  ('Flow Cytometry EC50', 'Primary Screening', '4PL logistic', 'Flow cytometry half-maximal effective concentration'),
  -- Epitope Binning
  ('Competitive SPR % Inhibition', 'Epitope Binning', 'Matrix thresholding', 'Percent inhibition by competitive SPR binning'),
  -- Affinity Characterization
  ('SPR ka', 'Affinity Characterization', '1:1 Langmuir', 'Association rate constant by surface plasmon resonance'),
  ('SPR kd', 'Affinity Characterization', '1:1 Langmuir', 'Dissociation rate constant by surface plasmon resonance'),
  ('SPR KD', 'Affinity Characterization', '1:1 Langmuir', 'Equilibrium dissociation constant by surface plasmon resonance'),
  -- Functional Assays
  ('Cell-based Reporter EC50', 'Functional Assays', '4PL logistic', 'EC50 from cell-based reporter gene assay'),
  ('Cell-based Reporter Emax', 'Functional Assays', '4PL logistic', 'Maximum efficacy from cell-based reporter gene assay'),
  -- Developability
  ('Thermal Stability Tm', 'Developability', 'Boltzmann sigmoid', 'Melting temperature by DSF or DSC'),
  ('Aggregation % Monomer', 'Developability', 'Gaussian peak fitting', 'Percent monomer by SEC-HPLC'),
  ('Aggregation Aggregate Size', 'Developability', 'Gaussian peak fitting', 'Aggregate hydrodynamic radius by DLS'),
  ('Viscosity', 'Developability', 'Nonlinear regression', 'Concentration-dependent viscosity measurement');

-- Grant permissions for new tables
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA biologics TO :LOGIN;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA biologics TO :LOGIN;
