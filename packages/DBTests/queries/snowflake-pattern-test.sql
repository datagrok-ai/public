--name: SnowflakePatternsAll
--connection: SnowflakeDBTests
--test: DBTests:expectTable(SnowflakePatternsAll(), OpenFile('System:AppData/DBTests/snowflake/data1-30.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data;
--end

--name: SnowflakeIntTypePatternNone
--connection: SnowflakeDBTests
--input: int id = 20
--test: DBTests:expectTable(SnowflakeIntTypePatternNone(), OpenFile('System:AppData/DBTests/snowflake/data20.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE id = @id;
--end

--name: SnowflakeStringTypeIntPatternOpMore
--connection: SnowflakeDBTests
--input: string id = ">28" {pattern: int}
--test: DBTests:expectTable(SnowflakeStringTypeIntPatternOpMore(), OpenFile('System:AppData/DBTests/snowflake/data29-30.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpMoreEq
--connection: SnowflakeDBTests
--input: string id = ">=29" {pattern: int}
--test: DBTests:expectTable(SnowflakeStringTypeIntPatternOpMoreEq(), OpenFile('System:AppData/DBTests/snowflake/data29-30.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpLessEq
--connection: SnowflakeDBTests
--input: string id = "<=1" {pattern: int}
--test: DBTests:expectTable(SnowflakeStringTypeIntPatternOpLessEq(), OpenFile('System:AppData/DBTests/snowflake/data1.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpLess
--connection: SnowflakeDBTests
--input: string id = "<2" {pattern: int}
--test: DBTests:expectTable(SnowflakeStringTypeIntPatternOpLess(), OpenFile('System:AppData/DBTests/snowflake/data1.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpIn
--connection: SnowflakeDBTests
--input: string id = "in(29, 30)" {pattern: int}
--test: DBTests:expectTable(SnowflakeStringTypeIntPatternOpIn(), OpenFile('System:AppData/DBTests/snowflake/data29-30.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpNotIn
--connection: SnowflakeDBTests
--input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
--test: DBTests:expectTable(SnowflakeStringTypeIntPatternOpNotIn(), OpenFile('System:AppData/DBTests/snowflake/data1-20.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpMinMax
--connection: SnowflakeDBTests
--input: string id = "min-max 29-30" {pattern: int}
--test: DBTests:expectTable(SnowflakeStringTypeIntPatternOpMinMax(), OpenFile('System:AppData/DBTests/snowflake/data29-30.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpNotEq
--connection: SnowflakeDBTests
--input: string id = "!=1" {pattern: int}
--test: DBTests:expectTable(SnowflakeStringTypeIntPatternOpNotEq(), OpenFile('System:AppData/DBTests/snowflake/data2-30.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @id(id)
--end

--name: SnowflakeDoubleTypePatternNone
--connection: SnowflakeDBTests
--input: double some_number = 510.32
--test: DBTests:expectTable(SnowflakeDoubleTypePatternNone(), OpenFile('System:AppData/DBTests/snowflake/data1.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE some_number = @some_number;
--end

--name: SnowflakeStringTypePatternDoubleOpMore
--connection: SnowflakeDBTests
--input: string some_number = ">975" {pattern: double}
--test: DBTests:expectTable(SnowflakeStringTypePatternDoubleOpMore(), OpenFile('System:AppData/DBTests/snowflake/data10,26.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @some_number(some_number);
--end

--name: SnowflakeStringTypePatternDoubleOpMoreEq
--connection: SnowflakeDBTests
--input: string some_number = ">=975" {pattern: double}
--test: DBTests:expectTable(SnowflakeStringTypePatternDoubleOpMoreEq(), OpenFile('System:AppData/DBTests/snowflake/data10,26.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @some_number(some_number);
--end

--name: SnowflakeStringTypePatternDoubleOpLess
--connection: SnowflakeDBTests
--input: string some_number = "<20" {pattern: double}
--test: DBTests:expectTable(SnowflakeStringTypePatternDoubleOpLess(), OpenFile('System:AppData/DBTests/snowflake/data5.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @some_number(some_number);
--end

--name: SnowflakeStringTypePatternDoubleOpLessEq
--connection: SnowflakeDBTests
--input: string some_number = "<=20" {pattern: double}
--test: DBTests:expectTable(SnowflakeStringTypePatternDoubleOpLessEq(), OpenFile('System:AppData/DBTests/snowflake/data5.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @some_number(some_number);
--end

--name: SnowflakeStringTypePatternStringOpContains
--connection: SnowflakeDBTests
--input: string first_name = "contains Z" {pattern: string}
--test: DBTests:expectTable(SnowflakeStringTypePatternStringOpContains(), OpenFile('System:AppData/DBTests/snowflake/data25.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @first_name(first_name);
--end

--name: SnowflakeStringTypePatternStringOpStartsWith
--connection: SnowflakeDBTests
--input: string first_name = "starts with W" {pattern: string}
--test: DBTests:expectTable(SnowflakeStringTypePatternStringOpStartsWith(), OpenFile('System:AppData/DBTests/snowflake/data23.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @first_name(first_name);
--end

--name: SnowflakeStringTypePatternStringOpEndsWith
--connection: SnowflakeDBTests
--input: string first_name = "ends with y" {pattern: string}
--test: DBTests:expectTable(SnowflakeStringTypePatternStringOpEndsWith(), OpenFile('System:AppData/DBTests/snowflake/data6,23,25.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @first_name(first_name);
--end

--name: SnowflakeStringTypePatternStringOpIn
--connection: SnowflakeDBTests
--input: string country = "in (Poland, Brazil)" {pattern: string}
--test: DBTests:expectTable(SnowflakeStringTypePatternStringOpIn(), OpenFile('System:AppData/DBTests/snowflake/data2,5,20.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @country(country);
--end

--name: SnowflakeStringTypePatternStringOpRegex
--connection: SnowflakeDBTests
--input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
--test: DBTests:expectTable(SnowflakeStringTypePatternStringOpRegex(), OpenFile('System:AppData/DBTests/snowflake/data9.d42')) //skip: GROK-12289
SELECT * FROM test.mock_data WHERE @email(email);
--end

--name: SnowflakePatternsAllParams
--connection: SnowflakeDBTests
--input: string first_name = "starts with p" {pattern: string}
--input: string id = ">1" {pattern :int}
--input: bool bool = false
--input: string email = "contains com" {pattern: string}
--input: string some_number = ">20" {pattern: double}
--input: string country = "in (Indonesia)" {pattern: string}
--input: string date = "before 1/1/2022" {pattern: datetime}
--test: DBTests:expectTable(SnowflakePatternsAllParams(), OpenFile("System:AppData/DBTests/snowflake/data13.d42")) //skip: GROK-12289
SELECT * FROM test.mock_data
WHERE @first_name(first_name)
  AND @id(id)
           AND bool = @bool
           AND @email(email)
           AND @some_number(some_number)
           AND @country(country)
           AND @date(date);
--end
