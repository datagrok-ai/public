CREATE TABLE binary_types (binary_data binary(5), varbinary_data varbinary(5));

INSERT INTO binary_types(binary_data, varbinary_data) TABLE (CAST( 'Datagrok' AS VARBINARY(5)),
    CAST(123456 AS VARBINARY(5)));

