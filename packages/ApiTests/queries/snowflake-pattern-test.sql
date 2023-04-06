--name: SnowflakePatternsAll
--connection: SnowflakeApiTests
--test: ApiTests:expectTable(SnowflakePatternsAll(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data1-30.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data;
--end

-- INT PATTERN

--name: SnowflakeIntTypePatternNone
--connection: SnowflakeApiTests
--input: int id = 20
--test: ApiTests:expectTable(SnowflakeIntTypePatternNone(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data20.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE id = @id;
--end

--name: SnowflakeStringTypeIntPatternOpMore
--connection: SnowflakeApiTests
--input: string id = ">28" {pattern: int}
--test: ApiTests:expectTable(SnowflakeStringTypeIntPatternOpMore(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data29-30.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpMoreEq
--connection: SnowflakeApiTests
--input: string id = ">=29" {pattern: int}
--test: ApiTests:expectTable(SnowflakeStringTypeIntPatternOpMoreEq(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data29-30.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpLessEq
--connection: SnowflakeApiTests
--input: string id = "<=1" {pattern: int}
--test: ApiTests:expectTable(SnowflakeStringTypeIntPatternOpLessEq(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data1.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpLess
--connection: SnowflakeApiTests
--input: string id = "<2" {pattern: int}
--test: ApiTests:expectTable(SnowflakeStringTypeIntPatternOpLess(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data1.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpIn
--connection: SnowflakeApiTests
--input: string id = "in(29, 30)" {pattern: int}
--test: ApiTests:expectTable(SnowflakeStringTypeIntPatternOpIn(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data29-30.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpNotIn
--connection: SnowflakeApiTests
--input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
--test: ApiTests:expectTable(SnowflakeStringTypeIntPatternOpNotIn(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data1-20.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpMinMax
--connection: SnowflakeApiTests
--input: string id = "min-max 29-30" {pattern: int}
--test: ApiTests:expectTable(SnowflakeStringTypeIntPatternOpMinMax(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data29-30.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpNotEq
--connection: SnowflakeApiTests
--input: string id = "!=1" {pattern: int}
--test: ApiTests:expectTable(SnowflakeStringTypeIntPatternOpNotEq(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data2-30.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @id(id)
--end


--DOUBLE PATTERN

--name: SnowflakeDoubleTypePatternNone
--connection: SnowflakeApiTests
--input: double some_number = 510.32
--test: ApiTests:expectTable(SnowflakeDoubleTypePatternNone(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data1.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE some_number = @some_number;
--end

--name: SnowflakeStringTypePatternDoubleOpMore
--connection: SnowflakeApiTests
--input: string some_number = ">975" {pattern: double}
--test: ApiTests:expectTable(SnowflakeStringTypePatternDoubleOpMore(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data10,26.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @some_number(some_number);
--end

--name: SnowflakeStringTypePatternDoubleOpMoreEq
--connection: SnowflakeApiTests
--input: string some_number = ">=975" {pattern: double}
--test: ApiTests:expectTable(SnowflakeStringTypePatternDoubleOpMoreEq(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data10,26.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @some_number(some_number);
--end

--name: SnowflakeStringTypePatternDoubleOpLess
--connection: SnowflakeApiTests
--input: string some_number = "<20" {pattern: double}
--test: ApiTests:expectTable(SnowflakeStringTypePatternDoubleOpLess(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data5.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @some_number(some_number);
--end

--name: SnowflakeStringTypePatternDoubleOpLessEq
--connection: SnowflakeApiTests
--input: string some_number = "<=20" {pattern: double}
--test: ApiTests:expectTable(SnowflakeStringTypePatternDoubleOpLessEq(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data5.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @some_number(some_number);
--end

--STRING PATTERN

--name: SnowflakeStringTypePatternStringOpContains
--connection: SnowflakeApiTests
--input: string first_name = "contains Z" {pattern: string}
--test: ApiTests:expectTable(SnowflakeStringTypePatternStringOpContains(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data25.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @first_name(first_name);
--end

--name: SnowflakeStringTypePatternStringOpStartsWith
--connection: SnowflakeApiTests
--input: string first_name = "starts with W" {pattern: string}
--test: ApiTests:expectTable(SnowflakeStringTypePatternStringOpStartsWith(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data23.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @first_name(first_name);
--end

--name: SnowflakeStringTypePatternStringOpEndsWith
--connection: SnowflakeApiTests
--input: string first_name = "ends with y" {pattern: string}
--test: ApiTests:expectTable(SnowflakeStringTypePatternStringOpEndsWith(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data6,23,25.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @first_name(first_name);
--end

--name: SnowflakeStringTypePatternStringOpIn
--connection: SnowflakeApiTests
--input: string country = "in (Poland, Brazil)" {pattern: string}
--test: ApiTests:expectTable(SnowflakeStringTypePatternStringOpIn(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data2,5,20.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @country(country);
--end

--name: SnowflakeStringTypePatternStringOpRegex
--connection: SnowflakeApiTests
--input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
--test: ApiTests:expectTable(SnowflakeStringTypePatternStringOpRegex(), OpenFile('System:AppData/ApiTests/datasets/tests/snowflake/data9.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @email(email);
--end

--name: SnowflakePatternsAllParams
--connection: SnowflakeApiTests
--input: string first_name = "starts with p" {pattern: string}
--input: string id = ">1" {pattern :int}
--input: bool bool = false
--input: string email = "contains com" {pattern: string}
--input: string some_number = ">20" {pattern: double}
--input: string country = "in (Indonesia)" {pattern: string}
--input: string date = "before 1/1/2022" {pattern: datetime}
--test: ApiTests:expectTable(SnowflakePatternsAllParams(), OpenFile("System:AppData/ApiTests/datasets/tests/snowflake/data13.d42")) //skip: GROK-12289
SELECT * FROM test.mock_data
WHERE @first_name(first_name)
  AND @id(id)
           AND bool = @bool
           AND @email(email)
           AND @some_number(some_number)
           AND @country(country)
           AND @date(date);
--end
