CREATE TABLE BYTEA_DATA (data BYTEA);

INSERT INTO BYTEA_DATA VALUES (pg_read_binary_file('/etc/file.txt'));