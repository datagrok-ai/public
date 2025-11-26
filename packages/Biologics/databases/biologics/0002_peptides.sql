CREATE TABLE biologics.peptides (
  id SERIAL PRIMARY KEY,
  identifier TEXT GENERATED ALWAYS AS ('GROKPEP-' || lpad(id::text, 6, '0')) STORED UNIQUE,
  name TEXT,
  helm TEXT,
  molecule_structure TEXT
);

ALTER TABLE biologics.assay_results
  ADD COLUMN IF NOT EXISTS peptide_id INTEGER REFERENCES biologics.peptides(id) ON UPDATE CASCADE;

CREATE INDEX IF NOT EXISTS biologics_assay_results_peptide_id_idx ON biologics.assay_results (peptide_id);