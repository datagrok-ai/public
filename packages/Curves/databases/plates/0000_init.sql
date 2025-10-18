CREATE SCHEMA IF NOT EXISTS plates;

CREATE TABLE plates.templates (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,
    description TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    created_by UUID
);

CREATE TABLE plates.properties (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,
    type TEXT NOT NULL CHECK (type IN ('int', 'double', 'bool', 'datetime', 'string', 'jsonb')),
    scope TEXT NOT NULL CHECK (scope IN ('plate', 'well')),
    nullable BOOLEAN NOT NULL DEFAULT TRUE,
    semType TEXT,
    userEditable BOOLEAN DEFAULT TRUE,
    units TEXT,
    min REAL,
    max REAL,
    choices TEXT,
    friendlyName TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE plates.template_plate_properties (
    template_id INTEGER NOT NULL REFERENCES plates.templates(id) ON DELETE CASCADE,
    property_id INTEGER NOT NULL REFERENCES plates.properties(id) ON DELETE CASCADE,
    is_required BOOLEAN NOT NULL DEFAULT FALSE,
    default_value TEXT,
    PRIMARY KEY (template_id, property_id)
);

-- Template well properties junction
CREATE TABLE plates.template_well_properties (
    template_id INTEGER NOT NULL REFERENCES plates.templates(id) ON DELETE CASCADE,
    property_id INTEGER NOT NULL REFERENCES plates.properties(id) ON DELETE CASCADE,
    is_required BOOLEAN NOT NULL DEFAULT FALSE,
    default_value TEXT,
    PRIMARY KEY (template_id, property_id)
);

-- Allowed values for properties
CREATE TABLE plates.property_allowed_values (
    id SERIAL PRIMARY KEY,
    property_id INTEGER NOT NULL REFERENCES plates.properties(id) ON DELETE CASCADE,
    value_string TEXT,
    value_float REAL,
    value_int INTEGER
);

-- Plate types
CREATE TABLE plates.plate_types (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,
    rows SMALLINT NOT NULL,
    cols SMALLINT NOT NULL,
    max_volume REAL
);

-- Plates
CREATE TABLE plates.plates (
    id SERIAL PRIMARY KEY,
    plate_type_id INTEGER NOT NULL REFERENCES plates.plate_types(id),
    barcode TEXT,
    description TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    created_by UUID
);

-- Analysis runs
CREATE TABLE plates.analysis_runs (
    id SERIAL PRIMARY KEY,
    plate_id INTEGER NOT NULL REFERENCES plates.plates(id) ON DELETE CASCADE,
    run_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    run_by UUID,
    analysis_type TEXT NOT NULL,
    groups TEXT[]
);

-- Plate details (plate-level properties)
CREATE TABLE plates.plate_details (
    plate_id INTEGER NOT NULL REFERENCES plates.plates(id),
    property_id INTEGER NOT NULL REFERENCES plates.properties(id),
    value_datetime TIMESTAMP WITH TIME ZONE,
    value_uuid uuid,
    value_num REAL,
    value_string TEXT,
    value_bool BOOLEAN,
    value_jsonb JSONB,
    analysis_run_id INTEGER REFERENCES plates.analysis_runs(id) ON DELETE SET NULL
);

-- Plate wells
CREATE TABLE plates.plate_wells (
    plate_id INTEGER NOT NULL REFERENCES plates.plates(id),
    row SMALLINT,
    col SMALLINT,
    details JSONB,
    PRIMARY KEY (plate_id, row, col)
);

-- Plate well values (well-level properties)
CREATE TABLE plates.plate_well_values (
    plate_id INTEGER NOT NULL REFERENCES plates.plates(id),
    row SMALLINT,
    col SMALLINT,
    property_id INTEGER REFERENCES plates.properties(id),
    value_num REAL,
    value_string TEXT,
    value_bool BOOLEAN,
    analysis_run_id INTEGER REFERENCES plates.analysis_runs(id) ON DELETE SET NULL
);

-- Analysis run parameters
CREATE TABLE plates.analysis_run_parameters (
    analysis_run_id INTEGER NOT NULL REFERENCES plates.analysis_runs(id) ON DELETE CASCADE,
    property_id INTEGER NOT NULL REFERENCES plates.properties(id),
    value_num REAL,
    value_string TEXT,
    value_bool BOOLEAN,
    value_datetime TIMESTAMP WITH TIME ZONE,
    value_jsonb JSONB,
    PRIMARY KEY (analysis_run_id, property_id)
);

CREATE TABLE plates.analysis_results (
    id SERIAL PRIMARY KEY,
    analysis_run_id INTEGER NOT NULL REFERENCES plates.analysis_runs(id) ON DELETE CASCADE,
    group_combination TEXT[],
    property_id INTEGER NOT NULL REFERENCES plates.properties(id),
    value_num REAL,
    value_string TEXT,
    value_bool BOOLEAN,
    value_datetime TIMESTAMP WITH TIME ZONE,
    value_jsonb JSONB 
);

-- Insert initial data
INSERT INTO plates.plate_types (id, name, rows, cols) VALUES 
    (1, 'Generic 96 wells', 8, 12),
    (2, 'Generic 384 wells', 16, 24),
    (3, 'Generic 1536 wells', 32, 48);

INSERT INTO plates.properties (id, name, type, units, scope) VALUES 
    (1000, 'Volume', 'double', 'uL', 'well'),
    (1001, 'Concentration', 'double', 'uM', 'well'),
    (1002, 'Sample', 'string', 'well'),
    (1003, 'Well Role', 'string', 'well');

INSERT INTO plates.property_allowed_values (property_id, value_string)
VALUES
    (1003, 'Empty'),
    (1003, 'Sample'),
    (1003, 'DMSO'),
    (1003, 'Low Control'),
    (1003, 'High Control');

-- Grant permissions
GRANT ALL PRIVILEGES ON SCHEMA plates TO CURRENT_USER;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA plates TO :LOGIN;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA plates TO :LOGIN;