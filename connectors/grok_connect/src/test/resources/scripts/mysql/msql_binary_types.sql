CREATE TABLE BINARY_TYPES(binary_type binary, varbinary_type varbinary, blob_type blob);

INSERT INTO BINARY_TYPES(binary_type, varbinary_type, blob_type) VALUES (BINARY('Hello'), BINARY('datagrok'),
                                                                         BINARY('Datagrok'));
