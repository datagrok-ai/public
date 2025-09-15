-- Step 1: Clear all existing data to ensure a clean slate.
-- This is the safest way to transition from the flawed old model.
TRUNCATE TABLE
    plates.semantic_types,
    plates.properties,
    plates.property_allowed_values,
    plates.plate_types,
    plates.layouts,
    plates.plates,
    plates.plate_details,
    plates.templates,
    plates.template_plate_properties,
    plates.template_well_properties,
    plates.plate_layout_wells,
    plates.plate_wells,
    plates.plate_well_values
RESTART IDENTITY CASCADE;


-- Step 2: Restructure the 'plates.properties' table.
-- This single command performs all necessary modifications.
ALTER TABLE plates.properties
    -- First, remove the old, incorrect global unique constraint on 'name'.
    DROP CONSTRAINT properties_name_key,

    -- Add the 'template_id' column to define ownership.
    ADD COLUMN template_id INTEGER REFERENCES plates.templates(id),

    -- Add the 'scope' column to explicitly define the context ('plate' or 'well').
    ADD COLUMN scope TEXT,

    -- Add a constraint to ensure the 'scope' column only contains valid values.
    ADD CONSTRAINT check_scope CHECK (scope IN ('plate', 'well')),

    -- Add the final, correct unique constraint. A property name must be
    -- unique for its scope WITHIN its parent template.
    ADD CONSTRAINT properties_template_scope_name_key UNIQUE (template_id, name, scope);


-- Step 3: Drop the now-redundant join tables.
-- The 'template_id' and 'scope' columns in 'properties' make these obsolete.
DROP TABLE IF EXISTS plates.template_plate_properties;
DROP TABLE IF EXISTS plates.template_well_properties;