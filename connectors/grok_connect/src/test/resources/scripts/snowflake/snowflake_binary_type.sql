CREATE TABLE binary_table (b BINARY);
INSERT INTO binary_table (b) SELECT TO_BINARY(HEX_ENCODE('Datagrok'), 'HEX');
