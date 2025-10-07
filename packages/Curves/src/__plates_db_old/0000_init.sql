-- Explains the meaning of a scalar property.
CREATE TABLE plates.semantic_types (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,    -- e.g., "Molecule", "Cell", "Tissue", "Organism", "Treatment", "Drug", "Image"
    description TEXT
);

-- Similar to core.properties or MolTrack.properties
-- Consider merging
CREATE TABLE plates.properties (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL, 
    nullable BOOLEAN NOT NULL DEFAULT TRUE,   -- used for validation of the submitted data
    type TEXT CHECK (type IN ('int', 'double', 'bool', 'datetime', 'string')),
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

-- Allowed values for properties that are strings, floats, or integers.
-- When allowed values are present, the property value must be one of the allowed values.
CREATE TABLE plates.property_allowed_values (
    id SERIAL PRIMARY KEY,
    property_id INTEGER NOT NULL REFERENCES plates.properties(id),
    value_string TEXT,
    value_float REAL,
    value_int INTEGER
);

-- maximum volume per well, material, shape (round bottom, conical).
CREATE TABLE plates.plate_types (
    id SERIAL PRIMARY KEY,
    name TEXT,
    rows SMALLINT,
    cols SMALLINT,
    max_volume REAL  -- microliters
);

-- layouts are reused across assay runs
CREATE TABLE plates.layouts (
    id SERIAL PRIMARY KEY,
    rows SMALLINT,
    cols SMALLINT,
    description TEXT,
    details JSONB,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    created_by UUID
);

-- represents an assay plate
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

-- User-defined plate details
-- Stored in the Plate.meta in TypeScript
CREATE TABLE plates.plate_details (
    plate_id INTEGER NOT NULL REFERENCES plates.plates(id),
    property_id INTEGER NOT NULL REFERENCES plates.properties(id),

    value_datetime TIMESTAMP WITH TIME ZONE,
    value_uuid uuid,
    value_num REAL,
    value_string TEXT,
    value_bool BOOLEAN,

    PRIMARY KEY (plate_id, property_id)
);

-- A template specifies a set of properties applicable on plate and well levels.
-- See also template_plate_properties and template_well_properties.
CREATE TABLE plates.templates (
    id SERIAL PRIMARY KEY,
    plate_layout_id INTEGER REFERENCES plates.plates(id),
    name TEXT,
    description TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    created_by UUID
);

CREATE TABLE plates.template_plate_properties (
    template_id INTEGER NOT NULL references plates.templates(id),
    property_id INTEGER NOT NULL references plates.properties(id),
    PRIMARY KEY (template_id, property_id)
);

CREATE TABLE plates.template_well_properties (
    template_id INTEGER NOT NULL references plates.templates(id),
    property_id INTEGER NOT NULL references plates.properties(id),
    PRIMARY KEY (template_id, property_id)
);


CREATE TABLE plates.plate_layout_wells (
    layout_id INTEGER NOT NULL references plates.plates(id),
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
  value_bool BOOLEAN
);


INSERT INTO plates.semantic_types (name) values ('Molecule');
INSERT INTO plates.semantic_types (name) values ('Solvent');
INSERT INTO plates.semantic_types (name) values ('URL');
INSERT INTO plates.semantic_types (name) values ('Image');

INSERT INTO plates.properties (id, name, type, units) values (1000, 'Volume', 'double', 'uL');
INSERT INTO plates.properties (id, name, type, units) values (1001, 'Concentration', 'double', 'uM');
INSERT INTO plates.properties (id, name, type) values (1002, 'Sample', 'string');
INSERT INTO plates.properties (id, name, type) values (1003, 'Well Role', 'string');

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
