CREATE TABLE numeric_type (
    number_value NUMBER(6, 2),
    small_value NUMBER(9),
    float_value FLOAT(126),
    binary_float_value BINARY_FLOAT,
    binary_double_value BINARY_DOUBLE

);

INSERT INTO numeric_type(number_value, small_value, float_value, binary_float_value, binary_double_value)
VALUES (1.987, 1123, 1.2e-4, 1.17549E-38F, 2.22507485850720E-308);
INSERT INTO numeric_type(number_value, small_value, float_value, binary_float_value, binary_double_value)
VALUES (1.9, 1, 0.55, binary_float_infinity, 0.2222);
INSERT INTO numeric_type(number_value, small_value, float_value, binary_float_value, binary_double_value)
VALUES (13.0, 1244124, 12, -binary_float_infinity, 2.2222);