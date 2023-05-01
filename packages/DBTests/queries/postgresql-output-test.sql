--name: PostgresqlArrayType
--connection: PostgreSQLDBTests
--test: DBTests:expectTable(PostgresqlArrayType(), OpenFile('System:AppData/DBTests/postgresql/postgres_output_array.d42'))
SELECT * FROM sal_emp;
--end

--name: PostgresqlBitString1
--connection: PostgreSQLDBTests
--test: DBTests:expectTable(PostgresqlBitString1(), OpenFile('System:AppData/DBTests/postgresql/postgres_output_bitstring_1.d42'))
SELECT * FROM test1;
--end

--name: PostgresqlBitString2
--connection: PostgreSQLDBTests
--test: DBTests:expectTable(PostgresqlBitString2(), OpenFile('System:AppData/DBTests/postgresql/postgres_output_bitstring_2.d42'))
SELECT * FROM test2;
--end

--name: PostgresqlComposite
--connection: PostgreSQLDBTests
--test: DBTests:expectTable(PostgresqlComposite(), OpenFile('System:AppData/DBTests/postgresql/postgres_output_composite.d42'))
SELECT * FROM on_hand;
--end

--name: PostgresqlDates
--connection: PostgreSQLDBTests
--test: DBTests:expectTable(PostgresqlDates(), OpenFile('System:AppData/DBTests/postgresql/postgres_output_dates.d42'))
SELECT * FROM dates;
--end

--name: PostgresqlJSONB
--connection: PostgreSQLDBTests
--test: DBTests:expectTable(PostgresqlJSONB(), OpenFile('System:AppData/DBTests/postgresql/postgres_output_jsonb.d42'))
SELECT * FROM jsonb_data;
--end

--name: PostgresqlNumeric
--connection: PostgreSQLDBTests
--test: DBTests:expectTable(PostgresqlNumeric(), OpenFile('System:AppData/DBTests/postgresql/postgres_output_numeric.d42'))
SELECT * FROM numeric_data;
--end

--name: PostgresqlDouble
--connection: PostgreSQLDBTests
--test: DBTests:expectTable(PostgresqlDouble(), OpenFile('System:AppData/DBTests/postgresql/postgres_output_doubles.d42'))
SELECT * FROM doubles;
--end

--name: PostgresqlReal
--connection: PostgreSQLDBTests
--test: DBTests:expectTable(PostgresqlReal(), OpenFile('System:AppData/DBTests/postgresql/postgres_output_real.d42'))
SELECT * FROM reals;
--end

--name: PostgresqlBigInt
--connection: PostgreSQLDBTests
--test: DBTests:expectTable(PostgresqlBigInt(), OpenFile('System:AppData/DBTests/postgresql/postgres_output_bigint.d42'))
SELECT * FROM bigint_data;
--end

--name: PostgresqlNumericPrecision
--connection: PostgreSQLDBTests
--test: DBTests:expectTable(PostgresqlNumericPrecision(), OpenFile("System:AppData/DBTests/postgresql/postgres_output_numeric_precision.d42"))
SELECT * FROM numeric_data_precision;
--end

--name: PostgresqlNumericPrecisionScale
--connection: PostgreSQLDBTests
--test: DBTests:expectTable(PostgresqlNumericPrecisionScale(), OpenFile("System:AppData/DBTests/postgresql/postgres_output_numeric_precision_scale.d42"))
SELECT * FROM numeric_data_precision_scale;
--end

--name: PostgresqlSmallSerial
--connection: PostgreSQLDBTests
--test: DBTests:expectTable(PostgresqlSmallSerial(), OpenFile('System:AppData/DBTests/postgresql/postgres_output_serial.d42'))
SELECT * FROM cars_small;
--end

--name: PostgresqlSerial
--connection: PostgreSQLDBTests
--test: DBTests:expectTable(PostgresqlSerial(), OpenFile('System:AppData/DBTests/postgresql/postgres_output_serial.d42'))
SELECT * FROM cars;
--end

--name: PostgresqlBigSerial
--connection: PostgreSQLDBTests
--test: DBTests:expectTable(PostgresqlBigSerial(), OpenFile('System:AppData/DBTests/postgresql/postgres_output_bigserial.d42'))
SELECT * FROM cars_big;
--end

--name: PostgresqlUUID
--connection: PostgreSQLDBTests
--test: DBTests:expectTable(PostgresqlUUID(), OpenFile('System:AppData/DBTests/postgresql/postgres_output_uuid.d42'))
SELECT * FROM uuid_data;
--end

--name: PostgresqlXML
--connection: PostgreSQLDBTests
--test: DBTests:expectTable(PostgresqlXML(), OpenFile('System:AppData/DBTests/postgresql/postgres_output_xml.d42'))
SELECT * FROM xml_data;
--end
