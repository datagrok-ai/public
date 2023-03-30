CREATE TABLE ARRAY_TYPE(array_type Array(UInt32), array_array_type Array(Array(UInt32)));
INSERT INTO ARRAY_TYPE (array_type, array_array_type) VALUES ([1, 2, 3, 4], [[24], [421, 12, 4], []]);
