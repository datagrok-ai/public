DROP TABLE IF EXISTS plates.properties CASCADE;

CREATE TABLE plates.properties (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL UNIQUE, 
    nullable BOOLEAN NOT NULL DEFAULT TRUE,
    type TEXT CHECK (type IN ('int', 'double', 'bool', 'datetime', 'string')),
    semType TEXT,
    userEditable BOOLEAN,
    units TEXT,
    min REAL,
    max REAL,
    step REAL,
    showSlider REAL,
    showPlusMinus REAL,
    inputType TEXT,
    category TEXT,
    format TEXT,
    choices TEXT,
    validators TEXT,
    friendlyName TEXT,
    options JSONB,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);