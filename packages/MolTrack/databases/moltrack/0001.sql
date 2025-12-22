ALTER TABLE moltrack.properties
  DROP CONSTRAINT properties_value_type_check;

ALTER TABLE moltrack.properties
  ADD CONSTRAINT properties_value_type_check
    CHECK (value_type IN ('int', 'double', 'datetime', 'uuid', 'string', 'bool'));