CREATE TABLE numeric_types (
    bigint_data bigint, int_data int, smallint_data smallint, tinyint_data tinyint, bit_data bit,
    decimal_data decimal(5, 2), numeric_data numeric(10, 5)
    );

INSERT INTO numeric_types(bigint_data, int_data, smallint_data, tinyint_data, bit_data, decimal_data, numeric_data)
VALUES (9223372036854775807, 2147483647, 32767, 123, 1, 123.22, 12345.12000);

CREATE TABLE float_types (float_data1 float(24), float_data2 float(53), real_data real);
INSERT INTO float_types(float_data1, float_data2, real_data) VALUES (- 1.79E+308,   -2.23E-308, 124124.23555);
INSERT INTO float_types(float_data1, float_data2, real_data) VALUES (34636.34661,  2.23E-308, 0.0);