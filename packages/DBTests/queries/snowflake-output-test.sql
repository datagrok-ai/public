--name: SnowflakeBinaryType
--connection: SnowflakeDBTests
--test: Dbtests:expectTable(SnowflakeBinaryType(), OpenFile('System:AppData/Dbtests/snowflake/snowflake_output_binary.d42')) //cat: Snowflake
SELECT to_char(B) AS B FROM test.public.binary_table;
--end

--name: SnowflakeDateTypes
--connection: SnowflakeDBTests
--test: Dbtests:expectTable(SnowflakeDateTypes(), OpenFile('System:AppData/Dbtests/snowflake/snowflake_output_dates.d42')) //cat: Snowflake
SELECT * FROM test.public.dates_table;
--end

--name: SnowflakeGeoType
--connection: SnowflakeDBTests
--test: Dbtests:expectTable(SnowflakeGeoType(), OpenFile('System:AppData/Dbtests/snowflake/snowflake_output_geo.d42')) //cat: Snowflake
SELECT * FROM test.public.geospatial_table;
--end

--name: SnowflakeNumericTypes
--connection: SnowflakeDBTests
--test: Dbtests:expectTable(SnowflakeNumericTypes(), OpenFile('System:AppData/Dbtests/snowflake/snowflake_output_numeric.d42')) //cat: Snowflake
SELECT * FROM test.public.test_number;
--end

--name: SnowflakeSemiStructuredTypes
--connection: SnowflakeDBTests
--test: Dbtests:expectTable(SnowflakeSemiStructuredTypes(), OpenFile('System:AppData/Dbtests/snowflake/snowflake_output_semi_structured.d42')) //cat: Snowflake
SELECT * FROM test.public.demonstration1;
--end
