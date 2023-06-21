CREATE TABLE INTEGER_TYPES (smallint_type smallint, int_type integer, bigint_type bigint);
INSERT INTO INTEGER_TYPES(smallint_type, int_type, bigint_type) VALUES (0, 2147483647, 9223372036854775807);
INSERT INTO INTEGER_TYPES(smallint_type, int_type, bigint_type) VALUES (-32768, -2147483648, -9223372036854775808);

CREATE TABLE FLOAT_TYPES (decimal_type decimal(38, 37), real_type real, double_precision_type double precision);
INSERT INTO FLOAT_TYPES(decimal_type, real_type, double_precision_type) VALUES (5.23553, -2E-10, 3.14);
INSERT INTO FLOAT_TYPES(decimal_type, real_type, double_precision_type) VALUES (0.9999, 214144412.2,
                                                                                -9223372036854775808.0);
