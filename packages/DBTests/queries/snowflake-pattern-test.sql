--name: SnowflakePatternsAll
--connection: SnowflakeDBTests
--test: Dbtests:expectTable(SnowflakePatternsAll(), OpenFile('System:AppData/Dbtests/common/data1-30.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data;
--end

--name: SnowflakeIntTypePatternNone
--connection: SnowflakeDBTests
--input: int id = 20
--test: Dbtests:expectTable(SnowflakeIntTypePatternNone(20), OpenFile('System:AppData/Dbtests/common/data20.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE id = @id;
--end

--name: SnowflakeStringTypeIntPatternOpMore
--connection: SnowflakeDBTests
--input: string id = ">28" {pattern: int}
--test: Dbtests:expectTable(SnowflakeStringTypeIntPatternOpMore('>28'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpMoreEq
--connection: SnowflakeDBTests
--input: string id = ">=29" {pattern: int}
--test: Dbtests:expectTable(SnowflakeStringTypeIntPatternOpMoreEq('>=29'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpLessEq
--connection: SnowflakeDBTests
--input: string id = "<=1" {pattern: int}
--test: Dbtests:expectTable(SnowflakeStringTypeIntPatternOpLessEq('<=1'), OpenFile('System:AppData/Dbtests/common/data1.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpLess
--connection: SnowflakeDBTests
--input: string id = "<2" {pattern: int}
--test: Dbtests:expectTable(SnowflakeStringTypeIntPatternOpLess('<2'), OpenFile('System:AppData/Dbtests/common/data1.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpIn
--connection: SnowflakeDBTests
--input: string id = "in(29, 30)" {pattern: int}
--test: Dbtests:expectTable(SnowflakeStringTypeIntPatternOpIn(id='in(29, 30)'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpNotIn
--connection: SnowflakeDBTests
--input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
--test: Dbtests:expectTable(SnowflakeStringTypeIntPatternOpNotIn(id='not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)'), OpenFile('System:AppData/Dbtests/common/data1-20.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpMinMax
--connection: SnowflakeDBTests
--input: string id = "min-max 29-30" {pattern: int}
--test: Dbtests:expectTable(SnowflakeStringTypeIntPatternOpMinMax(id='min-max 29-30'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE @id(id)
--end

--name: SnowflakeStringTypeIntPatternOpNotEq
--connection: SnowflakeDBTests
--input: string id = "!=1" {pattern: int}
--test: Dbtests:expectTable(SnowflakeStringTypeIntPatternOpNotEq('!=1'), OpenFile('System:AppData/Dbtests/common/data2-30.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE @id(id)
--end

--name: SnowflakeDoubleTypePatternNone
--connection: SnowflakeDBTests
--input: double some_number = 510.32
--test: Dbtests:expectTable(SnowflakeDoubleTypePatternNone(510.32), OpenFile('System:AppData/Dbtests/common/data1.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE some_number = @some_number;
--end

--name: SnowflakeStringTypePatternDoubleOpMore
--connection: SnowflakeDBTests
--input: string some_number = ">975" {pattern: double}
--test: Dbtests:expectTable(SnowflakeStringTypePatternDoubleOpMore('>975'), OpenFile('System:AppData/Dbtests/common/data10,26.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE @some_number(some_number);
--end

--name: SnowflakeStringTypePatternDoubleOpMoreEq
--connection: SnowflakeDBTests
--input: string some_number = ">=975" {pattern: double}
--test: Dbtests:expectTable(SnowflakeStringTypePatternDoubleOpMoreEq('>=975'), OpenFile('System:AppData/Dbtests/common/data10,26.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE @some_number(some_number);
--end

--name: SnowflakeStringTypePatternDoubleOpLess
--connection: SnowflakeDBTests
--input: string some_number = "<20" {pattern: double}
--test: Dbtests:expectTable(SnowflakeStringTypePatternDoubleOpLess('<20'), OpenFile('System:AppData/Dbtests/common/data5.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE @some_number(some_number);
--end

--name: SnowflakeStringTypePatternDoubleOpLessEq
--connection: SnowflakeDBTests
--input: string some_number = "<=20" {pattern: double}
--test: Dbtests:expectTable(SnowflakeStringTypePatternDoubleOpLessEq('<=20'), OpenFile('System:AppData/Dbtests/common/data5.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE @some_number(some_number);
--end

--name: SnowflakeStringTypePatternStringOpContains
--connection: SnowflakeDBTests
--input: string first_name = "contains Z" {pattern: string}
--test: Dbtests:expectTable(SnowflakeStringTypePatternStringOpContains(first_name = 'contains Z'), OpenFile('System:AppData/Dbtests/common/data25.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE @first_name(first_name);
--end

--name: SnowflakeStringTypePatternStringOpStartsWith
--connection: SnowflakeDBTests
--input: string first_name = "starts with W" {pattern: string}
--test: Dbtests:expectTable(SnowflakeStringTypePatternStringOpStartsWith(first_name='starts with W'), OpenFile('System:AppData/Dbtests/common/data23.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE @first_name(first_name);
--end

--name: SnowflakeStringTypePatternStringOpEndsWith
--connection: SnowflakeDBTests
--input: string first_name = "ends with y" {pattern: string}
--test: Dbtests:expectTable(SnowflakeStringTypePatternStringOpEndsWith(first_name = 'ends with y'), OpenFile('System:AppData/Dbtests/common/data6,23,25.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE @first_name(first_name);
--end

--name: SnowflakeStringTypePatternStringOpIn
--connection: SnowflakeDBTests
--input: string country = "in (Poland, Brazil)" {pattern: string}
--test: Dbtests:expectTable(SnowflakeStringTypePatternStringOpIn(country='in (Poland, Brazil)'), OpenFile('System:AppData/Dbtests/common/data2,5,20.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE @country(country);
--end

--name: SnowflakeStringTypePatternStringOpRegex
--connection: SnowflakeDBTests
--input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
--test: Dbtests:expectTable(SnowflakeStringTypePatternStringOpRegex(email = 'regex ^([A-Za-z0-9_]+@google.com.au)$'), OpenFile('System:AppData/Dbtests/common/data9.d42')) //cat: Snowflake
SELECT * FROM test.public.mock_data WHERE @email(email);
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
--test: Dbtests:expectTable(SnowflakePatternsAllParams(first_name = "starts with p", id = ">1", false, email = "contains com", some_number = ">20", country = "in (Indonesia)", date = "before 1/1/2022"), OpenFile("System:AppData/Dbtests/common/data13.d42")) //cat: Snowflake
SELECT * FROM test.public.mock_data
WHERE @first_name(first_name)
  AND @id(id)
           AND bool = @bool
           AND @email(email)
           AND @some_number(some_number)
           AND @country(country)
           AND @date(date);
--end
