DROP SCHEMA if exists plates cascade;
CREATE SCHEMA plates;

CREATE TABLE plates.plates (
    id SERIAL PRIMARY KEY,
    rows SMALLINT,
    cols SMALLINT,
    barcode TEXT,
    description TEXT,
    details JSONB,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    created_by UUID references users(id)
);

-- layouts are reused across assay runs
CREATE TABLE plates.layouts (
    id SERIAL PRIMARY KEY,
    rows SMALLINT,
    cols SMALLINT,
    description TEXT,
    details JSONB,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    created_by UUID references users(id)
);

CREATE TABLE plates.well_roles (
    id SERIAL PRIMARY KEY,
    name TEXT
);

CREATE TABLE plates.plate_layout_wells (
    layout_id INTEGER NOT NULL references plates.plates(id),
    row SMALLINT,
    col SMALLINT,
    well_role_id INTEGER NOT NULL references plates.well_roles(id),
    compound TEXT,
    volume FLOAT8,         -- microliters
    concentration FLOAT8,  -- micromolars
    PRIMARY KEY (layout_id, row, col)
);

CREATE TABLE plates.plate_wells (
    plate_id INTEGER NOT NULL references plates.plates(id),
    row SMALLINT,
    col SMALLINT,
    compound TEXT,
    details JSONB,
    PRIMARY KEY (plate_id, row, col)
);

INSERT INTO plates.well_roles (name) values ('Empty');
INSERT INTO plates.well_roles (name) values ('SAMPLE');
INSERT INTO plates.well_roles (name) values ('DMSO');
INSERT INTO plates.well_roles (name) values ('Low Control');
INSERT INTO plates.well_roles (name) values ('High Control');

GRANT ALL PRIVILEGES ON SCHEMA plates TO CURRENT_USER;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA plates TO CURRENT_USER;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA plates TO CURRENT_USER;