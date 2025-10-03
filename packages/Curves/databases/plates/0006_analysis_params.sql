-- Migration: Refactor to group-based analysis results system


ALTER TABLE plates.analysis_runs
    ADD COLUMN analysis_type TEXT,
    ADD COLUMN groups TEXT[];

UPDATE plates.analysis_runs SET analysis_type = analysis_name;

ALTER TABLE plates.analysis_runs
    ALTER COLUMN analysis_type SET NOT NULL;

ALTER TABLE plates.analysis_runs
    DROP COLUMN IF EXISTS parameters,
    DROP COLUMN IF EXISTS analysis_name;

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
    value_jsonb JSONB  -- for curve data
);

ALTER TABLE plates.plate_details 
    DROP COLUMN IF EXISTS group_type,
    DROP COLUMN IF EXISTS group_key;

DROP INDEX IF EXISTS plate_details_analysis_unique;
DROP INDEX IF EXISTS idx_plate_details_analysis;
DROP INDEX IF EXISTS idx_plate_details_group;

CREATE INDEX idx_analysis_results_run_id ON plates.analysis_results(analysis_run_id);
CREATE INDEX idx_analysis_results_groups ON plates.analysis_results USING GIN (group_combination);
CREATE INDEX idx_analysis_results_property ON plates.analysis_results(property_id);
CREATE INDEX idx_analysis_run_parameters_run_id ON plates.analysis_run_parameters(analysis_run_id);

COMMENT ON TABLE plates.analysis_runs IS 'Records analysis execution with discovered groups for the plate';
COMMENT ON COLUMN plates.analysis_runs.analysis_type IS 'Standardized analysis type identifier (DRC, qPCR, dose-ratio, etc.)';
COMMENT ON COLUMN plates.analysis_runs.groups IS 'All group values discovered in this analysis run (compounds, concentrations, etc.)';
COMMENT ON TABLE plates.analysis_run_parameters IS 'Analysis parameters using the properties system (regression type, confidence level, etc.)';
COMMENT ON TABLE plates.analysis_results IS 'Analysis results with group combinations and property values';
COMMENT ON COLUMN plates.analysis_results.group_combination IS 'Subset of run groups this result applies to, alphabetically sorted';

GRANT ALL PRIVILEGES ON TABLE plates.analysis_run_parameters TO :LOGIN;
GRANT ALL PRIVILEGES ON TABLE plates.analysis_results TO :LOGIN;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA plates TO :LOGIN;
GRANT ALL PRIVILEGES ON TABLE plates.analysis_runs TO :LOGIN;

