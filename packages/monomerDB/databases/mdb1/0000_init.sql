-- Create r_groups table first (no dependencies)
CREATE TABLE mdb1.r_groups (
    id SERIAL PRIMARY KEY,
    label VARCHAR(100), -- Added missing column referenced in index
    description TEXT,
    cap_group_smiles TEXT NOT NULL, -- SMILES of the capping group (e.g., O, [H], Cl...)
    cap_group_name VARCHAR(100), -- Name of the capping group (e.g., H, OH, Cl...)
    
    created_at TIMESTAMP NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMP NOT NULL DEFAULT NOW()
);

CREATE INDEX idx_r_groups_label ON mdb1.r_groups(label);

-- Populate r_groups table with common R-groups
INSERT INTO mdb1.r_groups (description, cap_group_smiles, cap_group_name)
VALUES
('Hydrogen', '[H]', 'H'),
('Hydroxyl', 'O', 'OH'),
('Amino', 'N', 'NH2'),
('Fluorine', 'F', 'F'),
('Chlorine', 'Cl', 'Cl'),
('Bromine', 'Br', 'Br'),
('Iodine', 'I', 'I');

-- Create monomer_libraries table second
CREATE TABLE mdb1.monomer_libraries (
    id SERIAL PRIMARY KEY,
    name VARCHAR(255) NOT NULL UNIQUE,
    friendly_name TEXT,
    description TEXT,
    library_type VARCHAR(50), -- e.g., 'internal', 'public', 'external'
    source VARCHAR(255), -- Original source or provider
    version VARCHAR(50),
    
    created_at TIMESTAMP NOT NULL DEFAULT NOW(),
    created_by VARCHAR(255) NOT NULL,
    updated_at TIMESTAMP NOT NULL DEFAULT NOW(),
    updated_by VARCHAR(255) NOT NULL
);

CREATE INDEX idx_monomer_libraries_name ON mdb1.monomer_libraries(name);
CREATE INDEX idx_monomer_libraries_type ON mdb1.monomer_libraries(library_type);

-- Create monomers table third (depends on r_groups and monomer_libraries)
CREATE TABLE mdb1.monomers (
    id SERIAL PRIMARY KEY,
    symbol VARCHAR(255) NOT NULL UNIQUE,
    name VARCHAR(255) NOT NULL,
    monomer_type VARCHAR(50) NOT NULL,
    polymer_type VARCHAR(50),
    grok_identifier TEXT GENERATED ALWAYS AS ('GROKMONO-' || id::text) STORED UNIQUE,
    
    -- Chemical structure
    smiles TEXT NOT NULL,
    molblock TEXT NOT NULL,
    capped_smiles TEXT,
    inchi_key TEXT, -- Added missing column referenced in index
    
    -- R-group references (nullable foreign keys)
    r1_id INTEGER REFERENCES mdb1.r_groups(id),
    r2_id INTEGER REFERENCES mdb1.r_groups(id),
    r3_id INTEGER REFERENCES mdb1.r_groups(id),
    r4_id INTEGER REFERENCES mdb1.r_groups(id),
    r5_id INTEGER REFERENCES mdb1.r_groups(id),
    r6_id INTEGER REFERENCES mdb1.r_groups(id),
    r7_id INTEGER REFERENCES mdb1.r_groups(id),
    r8_id INTEGER REFERENCES mdb1.r_groups(id),
    
    -- Library organization
    library_id INTEGER NOT NULL REFERENCES mdb1.monomer_libraries(id),
    
    -- Metadata
    created_at TIMESTAMP NOT NULL DEFAULT NOW(),
    created_by VARCHAR(255) NOT NULL,
    updated_at TIMESTAMP NOT NULL DEFAULT NOW(),
    updated_by VARCHAR(255) NOT NULL,
    
    -- Natural compound information
    natural_analog VARCHAR(255),
    -- Colors TBD
    
    CONSTRAINT unique_symbol_per_library UNIQUE(symbol, library_id)
);

CREATE INDEX idx_monomers_symbol ON mdb1.monomers(symbol);
CREATE INDEX idx_monomers_library ON mdb1.monomers(library_id);
CREATE INDEX idx_monomers_type ON mdb1.monomers(monomer_type);
CREATE INDEX idx_monomers_inchi_key ON mdb1.monomers(inchi_key);
CREATE INDEX idx_monomers_updated_at ON mdb1.monomers(updated_at);

-- Create monomer_properties table last (depends on monomers)
CREATE TABLE mdb1.monomer_properties (
    id SERIAL PRIMARY KEY,
    monomer_id INTEGER NOT NULL REFERENCES mdb1.monomers(id) ON DELETE CASCADE,
    name VARCHAR(255) NOT NULL,
    value_string TEXT,
    value_number REAL,
    value_boolean BOOLEAN,
    property_type VARCHAR(50), -- 'string', 'number', 'boolean', 'json'
    
    created_at TIMESTAMP NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMP NOT NULL DEFAULT NOW(),
    
    CONSTRAINT unique_property_per_monomer UNIQUE(monomer_id, name)
);

CREATE INDEX idx_monomer_properties_monomer ON mdb1.monomer_properties(monomer_id);
CREATE INDEX idx_monomer_properties_key ON mdb1.monomer_properties(name);

-- Grant permissions
GRANT ALL PRIVILEGES ON SCHEMA mdb1 TO CURRENT_USER;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA mdb1 TO :LOGIN;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA mdb1 TO :LOGIN;