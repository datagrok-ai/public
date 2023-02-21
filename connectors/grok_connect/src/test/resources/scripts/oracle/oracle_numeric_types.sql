CREATE TABLE numeric_type (
    number_value NUMBER(6, 2),
    small_value NUMBER(38),
    float_value FLOAT(7),
    binary_float_value BINARY_FLOAT,
    binary_double_value BINARY_DOUBLE

);

INSERT INTO numeric_type(number_value, small_value, float_value, binary_float_value, binary_double_value)
VALUES (9999.99, 100, 0.33, 3.40282E+38F, 1.79769313486231E+308);
INSERT INTO numeric_type(number_value, small_value, float_value, binary_float_value, binary_double_value)
VALUES (-9999.99, 1123, 1.2e-4, 1.17549E-38F, 2.22507485850720E-308);
INSERT INTO numeric_type(number_value, small_value, float_value, binary_float_value, binary_double_value)
VALUES (9.99, 1, 0.55, binary_float_infinity, 0.2222);
INSERT INTO numeric_type(number_value, small_value, float_value, binary_float_value, binary_double_value)
VALUES (13.00, 1244124, 12, -binary_float_infinity, 2.2222);