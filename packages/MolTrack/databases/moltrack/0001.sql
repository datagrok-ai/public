ALTER TABLE moltrack.properties
ADD COLUMN  min float, -- minimum value for numeric properties
ADD COLUMN  max float, -- maximum value for numeric properties
ADD COLUMN  choices text, -- JSON-encoded list of choices. Applicable to string properties only
ADD COLUMN  validators text -- JSON-encoded list of validators. Applicable to string properties only