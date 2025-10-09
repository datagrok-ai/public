CREATE TABLE plates.semantic_types (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,
    description TEXT
);

CREATE TABLE plates.templates (
    id SERIAL PRIMARY KEY,
    name TEXT,
    description TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    created_by UUID
);

CREATE TABLE plates.properties (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,
    nullable BOOLEAN NOT NULL DEFAULT TRUE,
    type TEXT CHECK (type IN ('int', 'double', 'bool', 'datetime', 'string')),
    semType TEXT,
    userEditable BOOLEAN,
    units TEXT,
    min REAL,
    max REAL,
    step REAL,
    showSlider REAL,
    showPlusMinus REAL,
    inputType TEXT,
    category TEXT,
    format TEXT,
    choices TEXT,
    validators TEXT,
    friendlyName TEXT,
    options JSONB,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    template_id INTEGER REFERENCES plates.templates(id),
    scope TEXT,
    
    CONSTRAINT check_scope CHECK (scope IN ('plate', 'well'))
);
CREATE UNIQUE INDEX properties_template_unique ON plates.properties(template_id, scope, name) WHERE template_id IS NOT NULL;

CREATE TABLE plates.property_allowed_values (
    id SERIAL PRIMARY KEY,
    property_id INTEGER NOT NULL REFERENCES plates.properties(id),
    value_string TEXT,
    value_float REAL,
    value_int INTEGER
);

CREATE TABLE plates.plate_types (
    id SERIAL PRIMARY KEY,
    name TEXT,
    rows SMALLINT,
    cols SMALLINT,
    max_volume REAL
);

CREATE TABLE plates.layouts (
    id SERIAL PRIMARY KEY,
    rows SMALLINT,
    cols SMALLINT,
    description TEXT,
    details JSONB,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    created_by UUID
);

CREATE TABLE plates.plates (
    id SERIAL PRIMARY KEY,
    layout_id INTEGER references plates.layouts(id),
    plate_type_id INTEGER references plates.plate_types(id),
    barcode TEXT,
    description TEXT,
    details JSONB,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    created_by UUID
);


CREATE TABLE plates.analysis_runs (
    id SERIAL PRIMARY KEY,
    plate_id INTEGER NOT NULL REFERENCES plates.plates(id) ON DELETE CASCADE,
    run_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    run_by UUID,
    analysis_type TEXT NOT NULL,
    groups TEXT[]
);
CREATE INDEX idx_analysis_runs_plate_id ON plates.analysis_runs(plate_id);
COMMENT ON COLUMN plates.analysis_runs.analysis_type IS 'Standardized analysis type identifier (DRC, qPCR, dose-ratio, etc.)';
COMMENT ON COLUMN plates.analysis_runs.groups IS 'All group values discovered in this analysis run (compounds, concentrations, etc.)';

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
CREATE UNIQUE INDEX plate_details_template_unique ON plates.plate_details(plate_id, property_id) WHERE analysis_run_id IS NULL;
CREATE INDEX idx_plate_details_jsonb ON plates.plate_details USING GIN (value_jsonb);
COMMENT ON COLUMN plates.plate_details.analysis_run_id IS 'Links to specific analysis run; NULL for template/manual data';

CREATE TABLE plates.plate_layout_wells (
    layout_id INTEGER NOT NULL references plates.layouts(id),
    row SMALLINT,
    col SMALLINT,
    PRIMARY KEY (layout_id, row, col)
);

CREATE TABLE plates.plate_wells (
    plate_id INTEGER NOT NULL references plates.plates(id),
    row SMALLINT,
    col SMALLINT,
    details JSONB,
    PRIMARY KEY (plate_id, row, col)
);

CREATE TABLE plates.plate_well_values (
    plate_id INTEGER NOT NULL references plates.plates(id),
    row SMALLINT,
    col SMALLINT,
    property_id INTEGER references plates.properties(id),
    value_num REAL,
    value_string TEXT,
    value_bool BOOLEAN,
    analysis_run_id INTEGER REFERENCES plates.analysis_runs(id) ON DELETE SET NULL
);
CREATE UNIQUE INDEX plate_well_values_template_unique ON plates.plate_well_values(plate_id, row, col, property_id) WHERE analysis_run_id IS NULL;
CREATE UNIQUE INDEX plate_well_values_analysis_unique ON plates.plate_well_values(plate_id, row, col, property_id, analysis_run_id) WHERE analysis_run_id IS NOT NULL;
CREATE INDEX idx_plate_well_values_analysis ON plates.plate_well_values(analysis_run_id) WHERE analysis_run_id IS NOT NULL;
COMMENT ON COLUMN plates.plate_well_values.analysis_run_id IS 'Links to specific analysis run; NULL for template/manual data';


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
CREATE INDEX idx_analysis_run_parameters_run_id ON plates.analysis_run_parameters(analysis_run_id);

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
CREATE INDEX idx_analysis_results_run_id ON plates.analysis_results(analysis_run_id);
CREATE INDEX idx_analysis_results_groups ON plates.analysis_results USING GIN (group_combination);
CREATE INDEX idx_analysis_results_property ON plates.analysis_results(property_id);
COMMENT ON COLUMN plates.analysis_results.group_combination IS 'Subset of run groups this result applies to, alphabetically sorted';


INSERT INTO plates.semantic_types (name) values ('Molecule'), ('Solvent'), ('URL'), ('Image');

INSERT INTO plates.properties (id, name, type, units, scope) values (1000, 'Volume', 'double', 'uL', 'well');
INSERT INTO plates.properties (id, name, type, units, scope) values (1001, 'Concentration', 'double', 'uM', 'well');
INSERT INTO plates.properties (id, name, type, scope) values (1002, 'Sample', 'string', 'well');
INSERT INTO plates.properties (id, name, type, scope) values (1003, 'Well Role', 'string', 'well');

INSERT INTO plates.plate_types (id, name, rows, cols) values (1, 'Generic 96 wells', 8, 12);
INSERT INTO plates.plate_types (id, name, rows, cols) values (2, 'Generic 384 wells', 16, 24);
INSERT INTO plates.plate_types (id, name, rows, cols) values (3, 'Generic 1536 wells', 32, 48);

INSERT INTO plates.property_allowed_values (property_id, value_string)
VALUES
    (1003, 'Empty'),
    (1003, 'Sample'),
    (1003, 'DMSO'),
    (1003, 'Low Control'),
    (1003, 'High Control');


GRANT ALL PRIVILEGES ON SCHEMA plates TO CURRENT_USER;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA plates TO :LOGIN;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA plates TO :LOGIN;