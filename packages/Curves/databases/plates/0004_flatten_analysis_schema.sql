-- Migration: Refactor to support analysis results in existing tables

DROP TABLE IF EXISTS plates.analysis_results_curves CASCADE;

ALTER TABLE plates.plate_details 
DROP CONSTRAINT plate_details_pkey,
ADD COLUMN group_type TEXT,
ADD COLUMN group_key TEXT;

CREATE UNIQUE INDEX plate_details_template_unique 
ON plates.plate_details(plate_id, property_id) 
WHERE analysis_run_id IS NULL;

CREATE UNIQUE INDEX plate_details_analysis_unique 
ON plates.plate_details(plate_id, property_id, analysis_run_id, group_key) 
WHERE analysis_run_id IS NOT NULL;

CREATE UNIQUE INDEX plate_well_values_template_unique 
ON plates.plate_well_values(plate_id, row, col, property_id) 
WHERE analysis_run_id IS NULL;

CREATE UNIQUE INDEX plate_well_values_analysis_unique 
ON plates.plate_well_values(plate_id, row, col, property_id, analysis_run_id) 
WHERE analysis_run_id IS NOT NULL;

CREATE INDEX idx_plate_details_analysis ON plates.plate_details(analysis_run_id, group_type, group_key) 
WHERE analysis_run_id IS NOT NULL;

CREATE INDEX idx_plate_details_group ON plates.plate_details(group_type, group_key) 
WHERE group_type IS NOT NULL;

CREATE INDEX idx_plate_well_values_analysis ON plates.plate_well_values(analysis_run_id) 
WHERE analysis_run_id IS NOT NULL;

COMMENT ON COLUMN plates.plate_details.group_type IS 'Type of grouping for analysis results (e.g., compound, concentration, timepoint)';
COMMENT ON COLUMN plates.plate_details.group_key IS 'Identifier within the group (e.g., Compound-A, 1e-7M, 24h)';
COMMENT ON COLUMN plates.plate_details.analysis_run_id IS 'Links to specific analysis run; NULL for template/manual data';
COMMENT ON COLUMN plates.plate_well_values.analysis_run_id IS 'Links to specific analysis run; NULL for template/manual data';

ALTER TABLE plates.properties
DROP CONSTRAINT IF EXISTS properties_name_key,
DROP CONSTRAINT IF EXISTS properties_template_scope_name_key;

CREATE UNIQUE INDEX properties_template_unique 
ON plates.properties(template_id, scope, name) 
WHERE template_id IS NOT NULL;

GRANT ALL PRIVILEGES ON TABLE plates.plate_details TO :LOGIN;
GRANT ALL PRIVILEGES ON TABLE plates.plate_well_values TO :LOGIN;