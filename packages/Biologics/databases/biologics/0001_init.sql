
ALTER TABLE biologics.assay_results ADD COLUMN IF NOT EXISTS adc_id INTEGER REFERENCES biologics.adc(id);

CREATE INDEX ON biologics.assay_results(adc_id);