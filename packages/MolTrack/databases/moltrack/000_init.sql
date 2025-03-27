-- Compounds table - core table for chemical structures
CREATE TABLE moltrack.compounds (
    id SERIAL PRIMARY KEY,
    canonical_smiles TEXT NOT NULL,
    inchi TEXT NOT NULL,
    inchikey TEXT NOT NULL UNIQUE,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    is_archived BOOLEAN DEFAULT FALSE
);

-- Batches table - for tracking physical samples
CREATE TABLE moltrack.batches (
    id SERIAL PRIMARY KEY,
    compound_id INTEGER NOT NULL REFERENCES moltrack.compounds(id),
    batch_number TEXT NOT NULL,
    amount REAL,
    amount_unit TEXT,
    purity REAL,
    vendor TEXT,
    catalog_id TEXT,
    acquisition_date DATE,
    expiry_date DATE,
    storage_location TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    created_by INTEGER,
    notes TEXT
);

-- Properties table - for calculated and measured properties
CREATE TABLE moltrack.properties (
    id SERIAL PRIMARY KEY,
    compound_id INTEGER NOT NULL REFERENCES moltrack.compounds(id),
    property_name TEXT NOT NULL,
    property_value TEXT NOT NULL,
    property_type TEXT CHECK (property_type IN ('CALCULATED', 'MEASURED', 'PREDICTED')),
    property_unit TEXT,
    calculation_method TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    created_by INTEGER
);


GRANT ALL ON TABLE moltrack.compounds to :LOGIN;
GRANT ALL ON TABLE moltrack.batches to :LOGIN;
GRANT ALL ON TABLE moltrack.properties to :LOGIN;



