-- forgot to add this (was initially saving as raw string)
ALTER TABLE plates.plate_details
ADD COLUMN value_jsonb JSONB;
GRANT ALL PRIVILEGES ON TABLE plates.plate_details TO :LOGIN;
CREATE INDEX idx_plate_details_jsonb ON plates.plate_details USING GIN (value_jsonb);