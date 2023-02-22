--name: PostgresqlArrayType
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlArrayType(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_array.csv'))
SELECT * FROM sal_emp;
--end

--name: PostgresqlBitString1
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlBitString1(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_bitstring_1.csv'))
SELECT * FROM test1;
--end

--name: PostgresqlBitString2
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlBitString2(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_bitstring_2.csv'))
SELECT * FROM test2;
--end

--name: PostgresqlComposite
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlComposite(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_composite.csv'))
SELECT * FROM on_hand;
--end

--name: PostgresqlDates
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlDates(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_dates.csv'))
SELECT * FROM dates;
--end

--name: PostgresqlJSONB
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlJSONB(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_jsonb.csv'))
SELECT * FROM jsonb_data;
--end

--name: PostgresqlNumeric
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlNumeric(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_numeric.csv'))
SELECT * FROM numeric_data;
--end

--name: PostgresqlDouble
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlDouble(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_doubles.csv'))
SELECT * FROM doubles;
--end

--name: PostgresqlReal
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlReal(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_real.csv'))
SELECT * FROM reals;
--end

--name: PostgresqlBigInt
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlBigInt(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_bigint.csv'))
SELECT * FROM bigint_data;
--end

--name: PostgresqlNumericPrecision
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlNumericPrecision(), OpenFile("System:AppData/ApiTests/datasets/tests/postgres_output_numeric_precision.csv"))
SELECT * FROM numeric_data_precision;
--end

--name: PostgresqlNumericPrecisionScale
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlNumericPrecisionScale(), OpenFile("System:AppData/ApiTests/datasets/tests/postgres_output_numeric_precision_scale.csv"))
SELECT * FROM numeric_data_precision_scale;
--end

--name: PostgresqlSmallSerial
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlSmallSerial(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_serial.csv'))
SELECT * FROM cars_small;
--end

--name: PostgresqlSerial
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlSerial(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_serial.csv'))
SELECT * FROM cars;
--end

--name: PostgresqlBigSerial
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlBigSerial(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_serial.csv'))
SELECT * FROM cars_big;
--end

--name: PostgresqlUUID
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlUUID(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_uuid.csv'))
SELECT * FROM uuid_data;
--end

--name: PostgresqlXML
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlXML(), OpenFile('System:AppData/ApiTests/datasets/tests/postgres_output_xml.csv'))
SELECT * FROM xml_data;
--end
