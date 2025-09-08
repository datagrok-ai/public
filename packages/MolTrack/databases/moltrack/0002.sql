CREATE TABLE moltrack.validators (
  id serial PRIMARY KEY,
  created_at timestamp with time zone DEFAULT (CURRENT_TIMESTAMP) NOT NULL,
  updated_at timestamp with time zone DEFAULT (CURRENT_TIMESTAMP) NOT NULL,
  created_by uuid NOT NULL REFERENCES moltrack.users (id),
  updated_by uuid NOT NULL REFERENCES moltrack.users (id),
  name text NOT NULL UNIQUE, -- e.g., "is_email", "is_url", "is_uuid", "is_smiles", "is_inchi", "is_inchikey"
  description text,
  entity_type text check (entity_type in ('BATCH', 'COMPOUND', 'ASSAY', 'ASSAY_RUN', 'ASSAY_RESULT')) NOT NULL,
  expression text NOT NULL
);