--name: SnowflakeBinaryType
--connection: SnowflakeDBTests
--test: Dbtests:expectTable(SnowflakeBinaryType(), OpenFile('System:AppData/Dbtests/snowflake/snowflake_output_binary.d42')) //skip: GROK-12289
SELECT to_char(B) AS B FROM test.binary_table;
--end

--name: SnowflakeDateTypes
--connection: SnowflakeDBTests
--test: Dbtests:expectTable(SnowflakeDateTypes(), OpenFile('System:AppData/Dbtests/snowflake/snowflake_output_dates.d42')) //skip: GROK-12289
SELECT * FROM test.dates_table;
--end

--name: SnowflakeGeoType
--connection: SnowflakeDBTests
--test: Dbtests:expectTable(SnowflakeGeoType(), OpenFile('System:AppData/Dbtests/snowflake/snowflake_output_geo.d42')) //skip: GROK-12289
SELECT * FROM test.geospatial_table;
--end

--name: SnowflakeNumericTypes
--connection: SnowflakeDBTests
--test: Dbtests:expectTable(SnowflakeNumericTypes(), OpenFile('System:AppData/Dbtests/snowflake/snowflake_output_numeric.d42')) //skip: GROK-12289
SELECT * FROM test.test_number;
--end

--name: SnowflakeSemiStructuredTypes
--connection: SnowflakeDBTests
--test: Dbtests:expectTable(SnowflakeSemiStructuredTypes(), OpenFile('System:AppData/Dbtests/snowflake/snowflake_output_semi_structured.d42')) //skip: GROK-12289
SELECT * FROM test.demonstration1;
--end
