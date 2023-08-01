--name: PostgresqlArrayType
--connection: PostgreSQLDBTests
--test: Dbtests:expectTable(PostgresqlArrayType(), OpenFile('System:AppData/Dbtests/postgresql/postgres_output_array.d42')) // cat: Postgresql
SELECT * FROM sal_emp;
--end

--name: PostgresqlBitString1
--connection: PostgreSQLDBTests
--test: Dbtests:expectTable(PostgresqlBitString1(), OpenFile('System:AppData/Dbtests/postgresql/postgres_output_bitstring_1.d42')) // cat: Postgresql
SELECT * FROM test1;
--end

--name: PostgresqlBitString2
--connection: PostgreSQLDBTests
--test: Dbtests:expectTable(PostgresqlBitString2(), OpenFile('System:AppData/Dbtests/postgresql/postgres_output_bitstring_2.d42')) // cat: Postgresql
SELECT * FROM test2;
--end

--name: PostgresqlComposite
--connection: PostgreSQLDBTests
--test: Dbtests:expectTable(PostgresqlComposite(), OpenFile('System:AppData/Dbtests/postgresql/postgres_output_composite.d42')) // cat: Postgresql
SELECT * FROM on_hand;
--end

--name: PostgresqlDates
--connection: PostgreSQLDBTests
--test: Dbtests:expectTable(PostgresqlDates(), OpenFile('System:AppData/Dbtests/postgresql/postgres_output_dates.d42')) // cat: Postgresql
SELECT * FROM dates;
--end

--name: PostgresqlJSONB
--connection: PostgreSQLDBTests
--test: Dbtests:expectTable(PostgresqlJSONB(), OpenFile('System:AppData/Dbtests/postgresql/postgres_output_jsonb.d42')) // cat: Postgresql
SELECT * FROM jsonb_data;
--end

--name: PostgresqlNumeric
--connection: PostgreSQLDBTests
--test: Dbtests:expectTable(PostgresqlNumeric(), OpenFile('System:AppData/Dbtests/postgresql/postgres_output_numeric.d42')) // cat: Postgresql
SELECT * FROM numeric_data;
--end

--name: PostgresqlDouble
--connection: PostgreSQLDBTests
--test: Dbtests:expectTable(PostgresqlDouble(), OpenFile('System:AppData/Dbtests/postgresql/postgres_output_doubles.d42')) // cat: Postgresql
SELECT * FROM doubles;
--end

--name: PostgresqlReal
--connection: PostgreSQLDBTests
--test: Dbtests:expectTable(PostgresqlReal(), OpenFile('System:AppData/Dbtests/postgresql/postgres_output_real.d42')) // cat: Postgresql
SELECT * FROM reals;
--end

--name: PostgresqlBigInt
--connection: PostgreSQLDBTests
--test: Dbtests:expectTable(PostgresqlBigInt(), OpenFile('System:AppData/Dbtests/postgresql/postgres_output_bigint.d42')) // cat: Postgresql
SELECT * FROM bigint_data;
--end

--name: PostgresqlNumericPrecision
--connection: PostgreSQLDBTests
--test: Dbtests:expectTable(PostgresqlNumericPrecision(), OpenFile("System:AppData/Dbtests/postgresql/postgres_output_numeric_precision.d42")) // cat: Postgresql
SELECT * FROM numeric_data_precision;
--end

--name: PostgresqlNumericPrecisionScale
--connection: PostgreSQLDBTests
--test: Dbtests:expectTable(PostgresqlNumericPrecisionScale(), OpenFile("System:AppData/Dbtests/postgresql/postgres_output_numeric_precision_scale.d42")) // cat: Postgresql
SELECT * FROM numeric_data_precision_scale;
--end

--name: PostgresqlSmallSerial
--connection: PostgreSQLDBTests
--test: Dbtests:expectTable(PostgresqlSmallSerial(), OpenFile('System:AppData/Dbtests/postgresql/postgres_output_serial.d42')) // cat: Postgresql
SELECT * FROM cars_small;
--end

--name: PostgresqlSerial
--connection: PostgreSQLDBTests
--test: Dbtests:expectTable(PostgresqlSerial(), OpenFile('System:AppData/Dbtests/postgresql/postgres_output_serial.d42')) // cat: Postgresql
SELECT * FROM cars;
--end

--name: PostgresqlBigSerial
--connection: PostgreSQLDBTests
--test: Dbtests:expectTable(PostgresqlBigSerial(), OpenFile('System:AppData/Dbtests/postgresql/postgres_output_bigserial.d42')) // cat: Postgresql
SELECT * FROM cars_big;
--end

--name: PostgresqlUUID
--connection: PostgreSQLDBTests
--test: Dbtests:expectTable(PostgresqlUUID(), OpenFile('System:AppData/Dbtests/postgresql/postgres_output_uuid.d42')) // cat: Postgresql
SELECT * FROM uuid_data;
--end

--name: PostgresqlXML
--connection: PostgreSQLDBTests
--test: Dbtests:expectTable(PostgresqlXML(), OpenFile('System:AppData/Dbtests/postgresql/postgres_output_xml.d42')) // cat: Postgresql
SELECT * FROM xml_data;
--end
