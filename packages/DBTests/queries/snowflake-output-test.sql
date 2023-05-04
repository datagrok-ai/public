--name: SnowflakeBinaryType
--connection: SnowflakeDBTests
--test: DBTests:expectTable(SnowflakeBinaryType(), OpenFile('System:AppData/DBTests/snowflake/snowflake_output_binary.d42')) //skip: GROK-12289
SELECT to_char(B) AS B FROM test.binary_table;
--end

--name: SnowflakeDateTypes
--connection: SnowflakeDBTests
--test: DBTests:expectTable(SnowflakeDateTypes(), OpenFile('System:AppData/DBTests/snowflake/snowflake_output_dates.d42')) //skip: GROK-12289
SELECT * FROM test.dates_table;
--end

--name: SnowflakeGeoType
--connection: SnowflakeDBTests
--test: DBTests:expectTable(SnowflakeGeoType(), OpenFile('System:AppData/DBTests/snowflake/snowflake_output_geo.d42')) //skip: GROK-12289
SELECT * FROM test.geospatial_table;
--end

--name: SnowflakeNumericTypes
--connection: SnowflakeDBTests
--test: DBTests:expectTable(SnowflakeNumericTypes(), OpenFile('System:AppData/DBTests/snowflake/snowflake_output_numeric.d42')) //skip: GROK-12289
SELECT * FROM test.test_number;
--end

--name: SnowflakeSemiStructuredTypes
--connection: SnowflakeDBTests
--test: DBTests:expectTable(SnowflakeSemiStructuredTypes(), OpenFile('System:AppData/DBTests/snowflake/snowflake_output_semi_structured.d42')) //skip: GROK-12289
SELECT * FROM test.demonstration1;
--end
