

-- Explains the meaning of a scalar property.
CREATE TABLE plts.semantic_types (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,    -- e.g., "Molecule", "Cell", "Tissue", "Organism", "Treatment", "Drug", "Image"
    description TEXT
);

-- Similar to core.properties or MolTrack.properties
-- Consider merging
CREATE TABLE plts.properties (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL, 
    type TEXT CHECK (type IN ('int', 'double', 'bool', 'datetime', 'string')),
    scope TEXT NOT NULL CHECK (scope IN ('plate', 'well')), -- whether the property is associated with a given well or is plate-wide
    nullable BOOLEAN NOT NULL DEFAULT TRUE,   -- used for validation of the submitted data
    semType TEXT,         -- Semantic type
    userEditable BOOLEAN, -- Whether the property should be editable via the UI
    units TEXT,           -- Units of measurement
    min REAL,             -- Minimum value. Applicable to numerical properties only
    max REAL,             -- Maximum value. Applicable to numerical properties only
    step REAL,            -- Step to be used in a slider. Only applies to numerical properties
    showSlider REAL,      -- Whether a slider appears next to the number input. Applies to numerical columns only
    showPlusMinus REAL,   -- Whether a plus/minus clicker appears next to the number input. Applies to numerical columns only
    inputType TEXT,       -- Property input type (see DG.InputType)
    category TEXT,        -- Corresponding category on the context panel
    format TEXT,          -- Value format, such as '0.000'
    choices TEXT,         -- JSON-encoded list of choices. Applicable to string properties only
    validators TEXT,      -- JSON-encoded list of validators, such as '[">42"]'
    friendlyName TEXT,    -- Custom field friendly name shown in [PropertyGrid]
    options JSONB,        -- Additional options
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);


CREATE TABLE plts.templates (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,
    description TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    created_by UUID
);

CREATE TABLE plts.template_plate_properties (
    template_id INTEGER NOT NULL REFERENCES plts.templates(id) ON DELETE CASCADE,
    property_id INTEGER NOT NULL REFERENCES plts.properties(id) ON DELETE CASCADE,
    is_required BOOLEAN NOT NULL DEFAULT FALSE,
    default_value TEXT,
    PRIMARY KEY (template_id, property_id)
);

-- Template well properties junction
CREATE TABLE plts.template_well_properties (
    template_id INTEGER NOT NULL REFERENCES plts.templates(id) ON DELETE CASCADE,
    property_id INTEGER NOT NULL REFERENCES plts.properties(id) ON DELETE CASCADE,
    is_required BOOLEAN NOT NULL DEFAULT FALSE,
    default_value TEXT,
    PRIMARY KEY (template_id, property_id)
);

-- Allowed values for properties
CREATE TABLE plts.property_allowed_values (
    id SERIAL PRIMARY KEY,
    property_id INTEGER NOT NULL REFERENCES plts.properties(id) ON DELETE CASCADE,
    value_string TEXT,
    value_float REAL,
    value_int INTEGER
);

-- Plate types
-- maximum volume per well, material, shape (round bottom, conical).
CREATE TABLE plts.plate_types (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,
    rows SMALLINT NOT NULL,
    cols SMALLINT NOT NULL,
    max_volume REAL --microliters
);
 
-- Plates

CREATE TABLE plts.plates (
    id SERIAL PRIMARY KEY,
    plate_type_id INTEGER NOT NULL REFERENCES plts.plate_types(id),
    barcode TEXT,
    description TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    created_by UUID
);

-- Analysis runs
CREATE TABLE plts.analysis_runs (
    id SERIAL PRIMARY KEY,
    plate_id INTEGER NOT NULL REFERENCES plts.plates(id) ON DELETE CASCADE,
    run_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    run_by UUID,
    analysis_type TEXT NOT NULL,
    groups TEXT[]
);

-- Plate details (plate-level properties)
CREATE TABLE plts.plate_details (
    plate_id INTEGER NOT NULL REFERENCES plts.plates(id),
    property_id INTEGER NOT NULL REFERENCES plts.properties(id),

    value_datetime TIMESTAMP WITH TIME ZONE,
    value_uuid uuid,
    value_num REAL,
    value_string TEXT,
    value_bool BOOLEAN,
    value_jsonb JSONB,

    analysis_run_id INTEGER REFERENCES plts.analysis_runs(id) ON DELETE SET NULL
);

-- Plate wells
CREATE TABLE plts.plate_wells (
    plate_id INTEGER NOT NULL REFERENCES plts.plates(id),
    row SMALLINT,
    col SMALLINT,
    details JSONB,
    PRIMARY KEY (plate_id, row, col)
);

-- Plate well values (well-level properties)
CREATE TABLE plts.plate_well_values (
    plate_id INTEGER NOT NULL REFERENCES plts.plates(id),
    row SMALLINT,
    col SMALLINT,
    property_id INTEGER REFERENCES plts.properties(id),
    value_num REAL,
    value_string TEXT,
    value_bool BOOLEAN,
    analysis_run_id INTEGER REFERENCES plts.analysis_runs(id) ON DELETE SET NULL
);

-- Analysis run parameters
CREATE TABLE plts.analysis_run_parameters (
    analysis_run_id INTEGER NOT NULL REFERENCES plts.analysis_runs(id) ON DELETE CASCADE,
    property_id INTEGER NOT NULL REFERENCES plts.properties(id),
    value_num REAL,
    value_string TEXT,
    value_bool BOOLEAN,
    value_datetime TIMESTAMP WITH TIME ZONE,
    value_jsonb JSONB,
    PRIMARY KEY (analysis_run_id, property_id)
);

CREATE TABLE plts.analysis_results (
    id SERIAL PRIMARY KEY,
    analysis_run_id INTEGER NOT NULL REFERENCES plts.analysis_runs(id) ON DELETE CASCADE,
    group_combination TEXT[],
    property_id INTEGER NOT NULL REFERENCES plts.properties(id),
    value_num REAL,
    value_string TEXT,
    value_bool BOOLEAN,
    value_datetime TIMESTAMP WITH TIME ZONE,
    value_jsonb JSONB 
);


INSERT INTO plts.semantic_types (name) values ('Molecule');
INSERT INTO plts.semantic_types (name) values ('Solvent');
INSERT INTO plts.semantic_types (name) values ('URL');
INSERT INTO plts.semantic_types (name) values ('Image');

INSERT INTO plts.properties (id, name, type, units, scope) VALUES 
    (1000, 'Volume', 'double', 'uL', 'well'),
    (1001, 'Concentration', 'double', 'uM', 'well'),
    (1002, 'Sample', 'string', '', 'well'),
    (1003, 'Well Role', 'string', '', 'well');

-- Insert initial data
INSERT INTO plts.plate_types (id, name, rows, cols) VALUES 
    (1, 'Generic 96 wells', 8, 12),
    (2, 'Generic 384 wells', 16, 24),
    (3, 'Generic 1536 wells', 32, 48);

INSERT INTO plts.property_allowed_values (property_id, value_string)
VALUES
    (1003, 'Empty'),
    (1003, 'Sample'),
    (1003, 'DMSO'),
    (1003, 'Low Control'),
    (1003, 'High Control');

-- Grant permissions
GRANT ALL PRIVILEGES ON SCHEMA plts TO CURRENT_USER;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA plts TO :LOGIN;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA plts TO :LOGIN;