--name: PostgresqlArrayType
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlArrayType(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/postgres_output_array.d42'))
SELECT * FROM sal_emp;
--end

--name: PostgresqlBitString1
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlBitString1(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/postgres_output_bitstring_1.d42'))
SELECT * FROM test1;
--end

--name: PostgresqlBitString2
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlBitString2(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/postgres_output_bitstring_2.d42'))
SELECT * FROM test2;
--end

--name: PostgresqlComposite
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlComposite(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/postgres_output_composite.d42'))
SELECT * FROM on_hand;
--end

--name: PostgresqlDates
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlDates(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/postgres_output_dates.d42'))
SELECT * FROM dates;
--end

--name: PostgresqlJSONB
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlJSONB(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/postgres_output_jsonb.d42'))
SELECT * FROM jsonb_data;
--end

--name: PostgresqlNumeric
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlNumeric(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/postgres_output_numeric.d42'))
SELECT * FROM numeric_data;
--end

--name: PostgresqlDouble
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlDouble(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/postgres_output_doubles.d42'))
SELECT * FROM doubles;
--end

--name: PostgresqlReal
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlReal(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/postgres_output_real.d42'))
SELECT * FROM reals;
--end

--name: PostgresqlBigInt
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlBigInt(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/postgres_output_bigint.d42'))
SELECT * FROM bigint_data;
--end

--name: PostgresqlNumericPrecision
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlNumericPrecision(), OpenFile("System:AppData/ApiTests/datasets/tests/postgresql/postgres_output_numeric_precision.d42"))
SELECT * FROM numeric_data_precision;
--end

--name: PostgresqlNumericPrecisionScale
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlNumericPrecisionScale(), OpenFile("System:AppData/ApiTests/datasets/tests/postgresql/postgres_output_numeric_precision_scale.d42"))
SELECT * FROM numeric_data_precision_scale;
--end

--name: PostgresqlSmallSerial
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlSmallSerial(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/postgres_output_serial.d42'))
SELECT * FROM cars_small;
--end

--name: PostgresqlSerial
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlSerial(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/postgres_output_serial.d42'))
SELECT * FROM cars;
--end

--name: PostgresqlBigSerial
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlBigSerial(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/postgres_output_bigserial.d42'))
SELECT * FROM cars_big;
--end

--name: PostgresqlUUID
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlUUID(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/postgres_output_uuid.d42'))
SELECT * FROM uuid_data;
--end

--name: PostgresqlXML
--connection: PostgreSQLApiTests
--test: ApiTests:expectTable(PostgresqlXML(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/postgres_output_xml.d42'))
SELECT * FROM xml_data;
--end
