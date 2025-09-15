TRUNCATE TABLE
    plates.plates,
    plates.plate_details,
    plates.plate_well_values,
    plates.plate_wells,
    plates.plate_layout_wells,
    plates.templates,
    plates.properties,
    plates.property_allowed_values,
    plates.layouts,
    plates.plate_types,
    plates.semantic_types
RESTART IDENTITY CASCADE;

-- Re-seed the basic data
INSERT INTO plates.semantic_types (name) values ('Molecule');
INSERT INTO plates.semantic_types (name) values ('Solvent');
INSERT INTO plates.semantic_types (name) values ('URL');
INSERT INTO plates.semantic_types (name) values ('Image');

INSERT INTO plates.plate_types (id, name, rows, cols) values (1, 'Generic 96 wells', 8, 12);
INSERT INTO plates.plate_types (id, name, rows, cols) values (2, 'Generic 384 wells', 16, 24);
INSERT INTO plates.plate_types (id, name, rows, cols) values (3, 'Generic 1536 wells', 32, 48);
