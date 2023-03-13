CREATE TABLE INTEGER_TYPES (tinyint_type tinyint, smallint_type smallint, mediumint_type mediumint,
int_type int, bigint_type bigint);

INSERT INTO INTEGER_TYPES(tinyint_type, smallint_type, mediumint_type, int_type, bigint_type)
VALUES (12, 32000, 167772, 2147483647, 9223372036854775807);

INSERT INTO INTEGER_TYPES(tinyint_type, smallint_type, mediumint_type, int_type, bigint_type)
VALUES (0, 1212, -1000, -2147483648, -9223372036854775808);
