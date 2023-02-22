--name: PostgresqlArrayType
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlArrayType(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_array.csv'))
SELECT * FROM sal_emp;
--end

--name: PostgresqlBitString1
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlBitString1(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_bitstring_1.csv'))
SELECT * FROM test1;
--end

--name: PostgresqlBitString2
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlBitString2(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_bitstring_2.csv'))
SELECT * FROM test2;
--end

--name: PostgresqlComposite
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlComposite(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_composite.csv'))
SELECT * FROM on_hand;
--end

--name: PostgresqlDates
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlDates(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_dates.csv'))
SELECT * FROM dates;
--end

--name: PostgresqlJSONB
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlJSONB(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_jsonb.csv'))
SELECT * FROM jsonb_data;
--end

--name: PostgresqlNumeric
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlNumeric(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_numeric.csv'))
SELECT * FROM numeric_data;
--end

--name: PostgresqlDouble
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlDouble(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_doubles.csv'))
SELECT * FROM doubles;
--end

--name: PostgresqlReal
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlReal(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_real.csv'))
SELECT * FROM reals;
--end

--name: PostgresqlBigInt
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlBigInt(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_bigint.csv'))
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
--test: ApiTests:compareTables(PostgresqlSmallSerial(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_serial.csv'))
SELECT * FROM cars_small;
--end

--name: PostgresqlSerial
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlSerial(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_serial.csv'))
SELECT * FROM cars;
--end

--name: PostgresqlBigSerial
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlBigSerial(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_serial.csv'))
SELECT * FROM cars_big;
--end

--name: PostgresqlUUID
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlUUID(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_uuid.csv'))
SELECT * FROM uuid_data;
--end

--name: PostgresqlXML
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlXML(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_xml.csv'))
SELECT * FROM xml_data;
--end
