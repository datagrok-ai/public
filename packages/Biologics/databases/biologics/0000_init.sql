

-- 1. Sequences
CREATE TABLE biologics.sequences (
  id SERIAL PRIMARY KEY,
  identifier TEXT GENERATED ALWAYS AS ('GROKSEQ-' || lpad(id::text, 6, '0')) STORED UNIQUE,
  name TEXT,
  sequence_type TEXT NOT NULL CHECK (sequence_type IN ('PROTEIN', 'DNA', 'RNA')),
  sequence TEXT NOT NULL
);

-- 2. Drugs
CREATE TABLE biologics.drugs (
  id SERIAL PRIMARY KEY,
  identifier TEXT GENERATED ALWAYS AS ('GROKMOL-' || lpad(id::text, 6, '0')) STORED UNIQUE,
  name TEXT,
  smiles TEXT NOT NULL
);

-- 3. Linkers
CREATE TABLE biologics.linkers (
  id SERIAL PRIMARY KEY,
  identifier TEXT GENERATED ALWAYS AS ('GROKLINKER-' || lpad(id::text, 6, '0')) STORED UNIQUE,
  linker_type TEXT NOT NULL CHECK (linker_type IN ('SMALL', 'PROTEIN')),
  linker_molecule_smiles TEXT,
  linker_sequence TEXT,
  CHECK (
    (linker_type = 'SMALL' AND linker_molecule_smiles IS NOT NULL AND linker_sequence IS NULL) OR
    (linker_type = 'PROTEIN' AND linker_sequence IS NOT NULL AND linker_molecule_smiles IS NULL)
  )
);

-- 4. ADC
CREATE TABLE biologics.adc (
  id SERIAL PRIMARY KEY,
  identifier TEXT GENERATED ALWAYS AS ('GROKADC-' || lpad(id::text, 6, '0')) STORED UNIQUE,
  name TEXT,
  antibody_id INTEGER NOT NULL REFERENCES biologics.sequences(id) ON UPDATE CASCADE,
  linker_id INTEGER NOT NULL REFERENCES biologics.linkers(id) ON UPDATE CASCADE,
  drug_id INTEGER NOT NULL REFERENCES biologics.drugs(id) ON UPDATE CASCADE,
  glyph TEXT  -- UTF-8 string (e.g. base64-encoded PNG)
);

-- 5. Assay Types
CREATE TABLE biologics.assay_types (
  id SERIAL PRIMARY KEY,
  name TEXT NOT NULL,
  description TEXT
);

-- 7. Target Organisms (created before assay_results for FK)
CREATE TABLE biologics.target_organisms (
  id SERIAL PRIMARY KEY,
  identifier TEXT GENERATED ALWAYS AS ('GROKORG-' || lpad(id::text, 6, '0')) STORED UNIQUE,
  name TEXT NOT NULL
);

-- 6. Assay Results
CREATE TABLE biologics.assay_results (
  id SERIAL PRIMARY KEY,
  assay_id INTEGER NOT NULL REFERENCES biologics.assay_types(id) ON UPDATE CASCADE ON DELETE RESTRICT,
  target_organism_id INTEGER REFERENCES biologics.target_organisms(id) ON UPDATE CASCADE ON DELETE SET NULL,
  result_value REAL,
  units TEXT,
  measured_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

-- 8. Purification Batches (minimal structure; extend later as needed)
CREATE TABLE biologics.purification_batches (
  id SERIAL PRIMARY KEY,
  identifier TEXT GENERATED ALWAYS AS ('GROKPUR-' || lpad(id::text, 6, '0')) STORED UNIQUE,
  sequence_id INTEGER REFERENCES biologics.sequences(id) ON UPDATE CASCADE,
  name TEXT,
  notes TEXT,
  created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

-- 9. Expression Batches (minimal structure; extend later as needed)
CREATE TABLE biologics.expression_batches (
  id SERIAL PRIMARY KEY,
  identifier TEXT GENERATED ALWAYS AS ('GROKEXP-' || lpad(id::text, 6, '0')) STORED UNIQUE,
  sequence_id INTEGER REFERENCES biologics.sequences(id) ON UPDATE CASCADE,
  expression_system TEXT,
  yield_mg REAL,
  name TEXT,
  notes TEXT,
  created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

-- Sample data: target organisms
INSERT INTO biologics.target_organisms (name) VALUES
  ('E. coli'),
  ('S. cerevisiae'),
  ('CHO');

-- Sample data: assay types
INSERT INTO biologics.assay_types (name, description) VALUES
  ('Caspase activity', 'Apoptosis induction via caspase activation readout'),
  ('IC50', 'Cytotoxic potency (half-maximal inhibitory concentration)'),
  ('Internalization half-life (t½)', 'Time to reach 50% internalization of the construct'),
  ('Cell binding constant (KD / EC50)', 'Flow cytometry binding curve-derived apparent affinity'),
  ('Average Drug-to-Antibody Ratio (DAR)', 'Average number of drug molecules per antibody (MS / HIC / UV/Vis)'),
  ('Binding affinity (KD / KDapp / EC50 / IC50)', 'Affinity from SPR, BLI, or ELISA dose-response experiments'),
  ('Cmax', 'Maximum observed concentration in PK profile'),
  ('Tmax', 'Time at which Cmax is observed'),
  ('AUC', 'Area under the concentration–time curve (systemic exposure)');

-- Helpful indexes
CREATE INDEX ON biologics.assay_results (assay_id);
CREATE INDEX ON biologics.assay_results (target_organism_id);
CREATE INDEX ON biologics.adc (antibody_id);
CREATE INDEX ON biologics.adc (linker_id);
CREATE INDEX ON biologics.adc (drug_id);

GRANT ALL PRIVILEGES ON SCHEMA biologics TO CURRENT_USER;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA biologics TO :LOGIN;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA biologics TO :LOGIN;
