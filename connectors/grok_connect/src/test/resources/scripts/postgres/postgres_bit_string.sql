CREATE TABLE test1 (a BIT(3), b BIT VARYING(5));
INSERT INTO test1 VALUES (B'101', B'0011');
INSERT INTO test1 VALUES (B'001', B'101');

CREATE TABLE test2 (c BIT(1), d BIT VARYING(5));
INSERT INTO test2 VALUES (B'1', B'101');
INSERT INTO test2 VALUES (B'0', B'0');