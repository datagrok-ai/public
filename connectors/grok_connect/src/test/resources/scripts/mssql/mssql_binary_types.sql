CREATE TABLE binary_types (binary_data binary(8), varbinary_data varbinary(6));
INSERT INTO binary_types(binary_data, varbinary_data) VALUES (CAST('Datagrok' AS BINARY(8)),
    CAST(123456 AS VARBINARY(6)));

