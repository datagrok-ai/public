-- Migration: Replace single sequence column with separate heavy_chain and light_chain columns,
-- annotate sequence regions with the chain they target, and add an assay_curves table.

-- 1. Sequences: add heavy_chain and light_chain, drop the single sequence column.
ALTER TABLE biologics.sequences
  ADD COLUMN IF NOT EXISTS heavy_chain TEXT,
  ADD COLUMN IF NOT EXISTS light_chain TEXT;

ALTER TABLE biologics.sequences
  DROP COLUMN IF EXISTS sequence;

-- sequence_type was tied to the dropped column; antibody chains are always proteins, so allow NULL.
ALTER TABLE biologics.sequences
  ALTER COLUMN sequence_type DROP NOT NULL;

-- 2. Sequence regions now target a specific chain (HC or LC); positions are within that chain.
ALTER TABLE biologics.sequence_regions
  ADD COLUMN IF NOT EXISTS chain TEXT CHECK (chain IN ('HC', 'LC'));

CREATE INDEX IF NOT EXISTS biologics_sequence_regions_chain_idx ON biologics.sequence_regions (chain);

-- 3. Sequence properties are computed per chain (HC and LC differ in MW, pI, charge, etc.).
ALTER TABLE biologics.sequence_properties
  ADD COLUMN IF NOT EXISTS chain TEXT CHECK (chain IN ('HC', 'LC'));

CREATE INDEX IF NOT EXISTS biologics_sequence_properties_chain_idx ON biologics.sequence_properties (chain);

-- 4. Liabilities: position is chain-relative, so the chain must be recorded alongside it.
ALTER TABLE biologics.sequence_liabilities
  ADD COLUMN IF NOT EXISTS chain TEXT CHECK (chain IN ('HC', 'LC'));

CREATE INDEX IF NOT EXISTS biologics_sequence_liabilities_chain_idx ON biologics.sequence_liabilities (chain);

-- 5. Expression and purification batches are produced per chain in real workflows.
ALTER TABLE biologics.expression_batches
  ADD COLUMN IF NOT EXISTS chain TEXT CHECK (chain IN ('HC', 'LC'));

CREATE INDEX IF NOT EXISTS biologics_expression_batches_chain_idx ON biologics.expression_batches (chain);

ALTER TABLE biologics.purification_batches
  ADD COLUMN IF NOT EXISTS chain TEXT CHECK (chain IN ('HC', 'LC'));

CREATE INDEX IF NOT EXISTS biologics_purification_batches_chain_idx ON biologics.purification_batches (chain);

-- 6. Assay curves: dose-response (or other fit) curves attached to assay results.
CREATE TABLE biologics.assay_curves (
  id SERIAL PRIMARY KEY,
  identifier TEXT GENERATED ALWAYS AS ('GROKCRV-' || lpad(id::text, 6, '0')) STORED UNIQUE,
  assay_result_id INTEGER NOT NULL REFERENCES biologics.assay_results(id) ON UPDATE CASCADE ON DELETE CASCADE,
  curve TEXT NOT NULL,
  created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX ON biologics.assay_curves (assay_result_id);

GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA biologics TO :LOGIN;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA biologics TO :LOGIN;
