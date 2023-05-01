--name: AthenaPatternsAll
--connection: AthenaDBTests
--test: DBTests:expectTable(AthenaPatternsAll(), OpenFile('System:AppData/DBTests/athena/data1-30.d42')) //skip:GROK-12423
SELECT * FROM mock_data
--end

-- INT PATTERN

--name: AthenaIntTypePatternNone
--connection: AthenaDBTests
--input: int id = 20
--test: DBTests:expectTable(AthenaIntTypePatternNone(), OpenFile('System:AppData/DBTests/athena/data20.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE id = @id;
--end

--name: AthenaStringTypeIntPatternOpMore
--connection: AthenaDBTests
--input: string id = ">28" {pattern: int}
--test: DBTests:expectTable(AthenaStringTypeIntPatternOpMore(), OpenFile('System:AppData/DBTests/athena/data29-30.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE @id(id)
--end

--name: AthenaStringTypeIntPatternOpMoreEq
--connection: AthenaDBTests
--input: string id = ">=29" {pattern: int}
--test: DBTests:expectTable(AthenaStringTypeIntPatternOpMoreEq(), OpenFile('System:AppData/DBTests/athena/data29-30.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE @id(id)
--end

--name: AthenaStringTypeIntPatternOpLessEq
--connection: AthenaDBTests
--input: string id = "<=1" {pattern: int}
--test: DBTests:expectTable(AthenaStringTypeIntPatternOpLessEq(), OpenFile('System:AppData/DBTests/athena/data1.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE @id(id)
--end

--name: AthenaStringTypeIntPatternOpLess
--connection: AthenaDBTests
--input: string id = "<2" {pattern: int}
--test: DBTests:expectTable(AthenaStringTypeIntPatternOpLess(), OpenFile('System:AppData/DBTests/athena/data1.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE @id(id)
--end

--name: AthenaStringTypeIntPatternOpIn
--connection: AthenaDBTests
--input: string id = "in(29, 30)" {pattern: int}
--test: DBTests:expectTable(AthenaStringTypeIntPatternOpIn(), OpenFile('System:AppData/DBTests/athena/data29-30.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE @id(id)
--end

--name: AthenaStringTypeIntPatternOpNotIn
--connection: AthenaDBTests
--input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
--test: DBTests:expectTable(AthenaStringTypeIntPatternOpNotIn(), OpenFile('System:AppData/DBTests/athena/data1-20.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE @id(id)
--end

--name: AthenaStringTypeIntPatternOpMinMax
--connection: AthenaDBTests
--input: string id = "min-max 29-30" {pattern: int}
--test: DBTests:expectTable(AthenaStringTypeIntPatternOpMinMax(), OpenFile('System:AppData/DBTests/athena/data29-30.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE @id(id)
--end

--name: AthenaStringTypeIntPatternOpNotEq
--connection: AthenaDBTests
--input: string id = "!=1" {pattern: int}
--test: DBTests:expectTable(AthenaStringTypeIntPatternOpNotEq(), OpenFile('System:AppData/DBTests/athena/data2-30.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE @id(id)
--end


--DOUBLE PATTERN

--name: AthenaDoubleTypePatternNone
--connection: AthenaDBTests
--input: double some_number = 510.32
--test: DBTests:expectTable(AthenaDoubleTypePatternNone(), OpenFile('System:AppData/DBTests/athena/data1.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE some_number = @some_number;
--end

--name: AthenaStringTypePatternDoubleOpMore
--connection: AthenaDBTests
--input: string some_number = ">975" {pattern: double}
--test: DBTests:expectTable(AthenaStringTypePatternDoubleOpMore(), OpenFile('System:AppData/DBTests/athena/data10,26.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: AthenaStringTypePatternDoubleOpMoreEq
--connection: AthenaDBTests
--input: string some_number = ">=975" {pattern: double}
--test: DBTests:expectTable(AthenaStringTypePatternDoubleOpMoreEq(), OpenFile('System:AppData/DBTests/athena/data10,26.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: AthenaStringTypePatternDoubleOpLess
--connection: AthenaDBTests
--input: string some_number = "<20" {pattern: double}
--test: DBTests:expectTable(AthenaStringTypePatternDoubleOpLess(), OpenFile('System:AppData/DBTests/athena/data5.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: AthenaStringTypePatternDoubleOpLessEq
--connection: AthenaDBTests
--input: string some_number = "<=20" {pattern: double}
--test: DBTests:expectTable(AthenaStringTypePatternDoubleOpLessEq(), OpenFile('System:AppData/DBTests/athena/data5.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--STRING PATTERN

--name: AthenaStringTypePatternStringOpContains
--connection: AthenaDBTests
--input: string first_name = "contains Z" {pattern: string}
--test: DBTests:expectTable(AthenaStringTypePatternStringOpContains(first_name = 'contains Z'), OpenFile('System:AppData/DBTests/athena/data25.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: AthenaStringTypePatternStringOpStartsWith
--connection: AthenaDBTests
--input: string first_name = "starts with W" {pattern: string}
--test: DBTests:expectTable(AthenaStringTypePatternStringOpStartsWith(first_name='starts with W'), OpenFile('System:AppData/DBTests/athena/data23.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: AthenaStringTypePatternStringOpEndsWith
--connection: AthenaDBTests
--input: string first_name = "ends with y" {pattern: string}
--test: DBTests:expectTable(AthenaStringTypePatternStringOpEndsWith(first_name = 'ends with y'), OpenFile('System:AppData/DBTests/athena/data6,23,25.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: AthenaStringTypePatternStringOpIn
--connection: AthenaDBTests
--input: string country = "in (Poland, Brazil)" {pattern: string}
--test: DBTests:expectTable(AthenaStringTypePatternStringOpIn(), OpenFile('System:AppData/DBTests/athena/data2,5,20.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE @country(country);
--end

--name: AthenaStringTypePatternStringOpRegex
--connection: AthenaDBTests
--input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
--test: DBTests:expectTable(AthenaStringTypePatternStringOpRegex(email = 'regex ^([A-Za-z0-9_]+@google.com.au)$'), OpenFile('System:AppData/DBTests/athena/data9.d42')) //skip:GROK-12423
SELECT * FROM mock_data WHERE @email(email);
--end

--name: AthenaPatternsAllParams
--connection: AthenaDBTests
--input: string first_name = "starts with p" {pattern: string}
--input: string id = ">1" {pattern :int}
--input: bool bool = false
--input: string email = "contains com" {pattern: string}
--input: string some_number = ">20" {pattern: double}
--input: string country = "in (Indonesia)" {pattern: string}
--input: string date = "before 1/1/2022" {pattern: datetime}
--test: DBTests:expectTable(AthenaPatternsAllParams(first_name = "starts with p", id = ">1", false, email = "contains com", some_number = ">20", country = "in (Indonesia)", date = "before 1/1/2022"), OpenFile("System:AppData/DBTests/athena/data13.d42")) //skip:GROK-12423
SELECT * FROM mock_data
WHERE @first_name(first_name)
  AND @id(id)
           AND bool = @bool
           AND @email(email)
           AND @some_number(some_number)
           AND @country(country)
           AND @date(date);
--end
