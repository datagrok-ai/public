--name: PostgresqlArrayType
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlArrayType(), ApiTests:getTable('postgres_output_array.csv', path='ApiTests/datasets/tests'))
SELECT * FROM sal_emp;
--end

--name: PostgresqlBitString1
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlBitString1(), ApiTests:getTable('postgres_output_bitstring_1.csv', path='ApiTests/datasets/tests'))
SELECT * FROM test1;
--end

--name: PostgresqlBitString2
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlBitString2(), ApiTests:getTable('postgres_output_bitstring_2.csv', path='ApiTests/datasets/tests'))
SELECT * FROM test2;
--end

--name: PostgresqlComposite
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlComposite(), ApiTests:getTable('postgres_output_composite.csv', path='ApiTests/datasets/tests'))
SELECT * FROM on_hand;
--end

--name: PostgresqlDates
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlDates(), ApiTests:getTable('postgres_output_dates.csv', path='ApiTests/datasets/tests'))
SELECT * FROM dates;
--end

--name: PostgresqlJSONB
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlJSONB(), ApiTests:getTable('postgres_output_jsonb.csv', path='ApiTests/datasets/tests'))
SELECT * FROM jsonb_data;
--end

--name: PostgresqlNumeric
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlNumeric(), ApiTests:getTable('postgres_output_numeric.csv', path='ApiTests/datasets/tests'))
SELECT * FROM numeric_data;
--end

--name: PostgresqlDouble
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlDouble(), ApiTests:getTable('postgres_output_doubles.csv', path='ApiTests/datasets/tests'))
SELECT * FROM doubles;
--end

--name: PostgresqlReal
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlReal(), ApiTests:getTable('postgres_output_real.csv', path='ApiTests/datasets/tests'))
SELECT * FROM reals;
--end

--name: PostgresqlBigInt
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlBigInt(), ApiTests:getTable('postgres_output_bigint.csv', path='ApiTests/datasets/tests'))
SELECT * FROM bigint_data;
--end

--name: PostgresqlNumericPrecision
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlNumericPrecision(), ApiTests:getTable("postgres_output_numeric_precision.csv", path="ApiTests/datasets/tests"))
SELECT * FROM numeric_data_precision;
--end

--name: PostgresqlNumericPrecisionScale
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlNumericPrecisionScale(), ApiTests:getTable("postgres_output_numeric_precision_scale.csv", path="ApiTests/datasets/tests"))
SELECT * FROM numeric_data_precision_scale;
--end

--name: PostgresqlSmallSerial
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlSmallSerial(), ApiTests:getTable('postgres_output_serial.csv', path='ApiTests/datasets/tests'))
SELECT * FROM cars_small;
--end

--name: PostgresqlSerial
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlSerial(), ApiTests:getTable('postgres_output_serial.csv', path='ApiTests/datasets/tests'))
SELECT * FROM cars;
--end

--name: PostgresqlBigSerial
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlBigSerial(), ApiTests:getTable('postgres_output_serial.csv', path='ApiTests/datasets/tests'))
SELECT * FROM cars_big;
--end

--name: PostgresqlUUID
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlUUID(), ApiTests:getTable('postgres_output_uuid.csv', path='ApiTests/datasets/tests'))
SELECT * FROM uuid_data;
--end

--name: PostgresqlXML
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlXML(), ApiTests:getTable('postgres_output_xml.csv', path='ApiTests/datasets/tests'))
SELECT * FROM xml_data;
--end
