-- Migration: Add support for analysis runs, curve results, and property origins.

-- This table records each execution of an analysis on a specific plate.
CREATE TABLE plates.analysis_runs (
    id SERIAL PRIMARY KEY,
    plate_id INTEGER NOT NULL REFERENCES plates.plates(id) ON DELETE CASCADE,
    analysis_name TEXT NOT NULL,
    run_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    run_by UUID,
    parameters JSONB -- Stores input parameters for the run, like mappings.
);
CREATE INDEX idx_analysis_runs_plate_id ON plates.analysis_runs(plate_id);
COMMENT ON TABLE plates.analysis_runs IS 'Records an execution per specific analysis type per specific plate.';


-- This table stores individual curve results from an analysis run for granular querying.
CREATE TABLE plates.analysis_results_curves (
    id SERIAL PRIMARY KEY,
    analysis_run_id INTEGER NOT NULL REFERENCES plates.analysis_runs(id) ON DELETE CASCADE,
    series_name TEXT,     -- eg Sample A-123, Control 
    curve_json JSONB,  
    ic50 REAL,
    hill_slope REAL,
    r_squared REAL,
    min_value REAL,
    max_value REAL,
    auc REAL,
    details JSONB       
);
CREATE INDEX idx_results_curves_run_id ON plates.analysis_results_curves(analysis_run_id);
COMMENT ON TABLE plates.analysis_results_curves IS 'Stores individual curve results from an analysis run.';

---

-- Link other analysis results back to their specific run.
-- This is for future-proofing if you have analyses that produce simple plate/well-level scalars.
ALTER TABLE plates.plate_details ADD COLUMN analysis_run_id INTEGER REFERENCES plates.analysis_runs(id) ON DELETE SET NULL;
ALTER TABLE plates.plate_well_values ADD COLUMN analysis_run_id INTEGER REFERENCES plates.analysis_runs(id) ON DELETE SET NULL;

COMMENT ON COLUMN plates.plate_details.analysis_run_id IS 'FK to analysis_runs. Links plate-level results to the run that produced them.';
COMMENT ON COLUMN plates.plate_well_values.analysis_run_id IS 'FK to analysis_runs. Links well-level results to the run that produced them.';

---

-- Add the 'origin_plate_id' bookmark to the properties table.
-- For "global" properties, this stores the ID of the plate that first created this property.
ALTER TABLE plates.properties ADD COLUMN origin_plate_id INTEGER REFERENCES plates.plates(id) ON DELETE SET NULL;
COMMENT ON COLUMN plates.properties.origin_plate_id IS 'For global properties, stores the ID of the plate that first created this property.';

GRANT ALL PRIVILEGES ON TABLE plates.analysis_runs TO :LOGIN;
GRANT ALL PRIVILEGES ON TABLE plates.analysis_results_curves TO :LOGIN;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA plates TO :LOGIN;