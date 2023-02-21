CREATE OR REPLACE TYPE mem_type IS VARRAY(10) of VARCHAR2(15);

CREATE TABLE varrays (members mem_type);

INSERT INTO varrays(members) VALUES (mem_type('Brenda','Richard'));
