--name: SnowflakeBinaryType
--connection: SnowflakeApiTests
--test: ApiTests:expectTable(SnowflakeBinaryType(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/snowflake_output_binary.d42')) //skip: GROK-12289
SELECT to_char(B) AS B FROM binary_table;
--end

--name: SnowflakeDateTypes
--connection: SnowflakeApiTests
--test: ApiTests:expectTable(SnowflakeDateTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/snowflake_output_dates.d42')) //skip: GROK-12289
SELECT * FROM dates_table;
--end

--name: SnowflakeGeoType
--connection: SnowflakeApiTests
--test: ApiTests:expectTable(SnowflakeGeoType(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/snowflake_output_geo.d42')) //skip: GROK-12289
SELECT * FROM geospatial_table;
--end

--name: SnowflakeNumericTypes
--connection: SnowflakeApiTests
--test: ApiTests:expectTable(SnowflakeNumericTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/snowflake_output_numeric.d42')) //skip: GROK-12289
SELECT * FROM test_number;
--end

--name: SnowflakeSemiStructuredTypes
--connection: SnowflakeApiTests
--test: ApiTests:expectTable(SnowflakeSemiStructuredTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/snowflake_output_semi_structured.d42')) //skip: GROK-12289
SELECT * FROM demonstration1;
--end
