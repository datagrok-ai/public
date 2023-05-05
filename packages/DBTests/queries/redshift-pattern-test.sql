--name: RedshiftPatternsAll
--connection: RedshiftDBTests
--test: DBTests:expectTable(RedshiftPatternsAll(), OpenFile('System:AppData/DBTests/redshift/data1-30.d42')) //skip: GROK-12507
SELECT * FROM mock_data
--end

--name: RedshiftIntTypePatternNone
--connection: RedshiftDBTests
--input: int id = 20
--test: DBTests:expectTable(RedshiftIntTypePatternNone(20), OpenFile('System:AppData/DBTests/redshift/data20.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE id = @id;
--end

--name: RedshiftStringTypeIntPatternOpMore
--connection: RedshiftDBTests
--input: string id = ">28" {pattern: int}
--test: DBTests:expectTable(RedshiftStringTypeIntPatternOpMore(), OpenFile('System:AppData/DBTests/redshift/data29-30.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @id(id)
--end

--name: RedshiftStringTypeIntPatternOpMoreEq
--connection: RedshiftDBTests
--input: string id = ">=29" {pattern: int}
--test: DBTests:expectTable(RedshiftStringTypeIntPatternOpMoreEq(), OpenFile('System:AppData/DBTests/redshift/data29-30.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @id(id)
--end

--name: RedshiftStringTypeIntPatternOpLessEq
--connection: RedshiftDBTests
--input: string id = "<=1" {pattern: int}
--test: DBTests:expectTable(RedshiftStringTypeIntPatternOpLessEq(), OpenFile('System:AppData/DBTests/redshift/data1.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @id(id)
--end

--name: RedshiftStringTypeIntPatternOpLess
--connection: RedshiftDBTests
--input: string id = "<2" {pattern: int}
--test: DBTests:expectTable(RedshiftStringTypeIntPatternOpLess(), OpenFile('System:AppData/DBTests/redshift/data1.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @id(id)
--end

--name: RedshiftStringTypeIntPatternOpIn
--connection: RedshiftDBTests
--input: string id = "in(29, 30)" {pattern: int}
--test: DBTests:expectTable(RedshiftStringTypeIntPatternOpIn(), OpenFile('System:AppData/DBTests/redshift/data29-30.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @id(id)
--end

--name: RedshiftStringTypeIntPatternOpNotIn
--connection: RedshiftDBTests
--input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
--test: DBTests:expectTable(RedshiftStringTypeIntPatternOpNotIn(), OpenFile('System:AppData/DBTests/redshift/data1-20.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @id(id)
--end

--name: RedshiftStringTypeIntPatternOpMinMax
--connection: RedshiftDBTests
--input: string id = "min-max 29-30" {pattern: int}
--test: DBTests:expectTable(RedshiftStringTypeIntPatternOpMinMax(), OpenFile('System:AppData/DBTests/redshift/data29-30.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @id(id)
--end

--name: RedshiftStringTypeIntPatternOpNotEq
--connection: RedshiftDBTests
--input: string id = "!=1" {pattern: int}
--test: DBTests:expectTable(RedshiftStringTypeIntPatternOpNotEq(), OpenFile('System:AppData/DBTests/redshift/data2-30.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @id(id)
--end

--name: RedshiftDoubleTypePatternNone
--connection: RedshiftDBTests
--input: double some_number = 510.32
--test: DBTests:expectTable(RedshiftDoubleTypePatternNone(510.32), OpenFile('System:AppData/DBTests/redshift/data1.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE some_number = @some_number;
--end

--name: RedshiftStringTypePatternDoubleOpMore
--connection: RedshiftDBTests
--input: string some_number = ">975" {pattern: double}
--test: DBTests:expectTable(RedshiftStringTypePatternDoubleOpMore(), OpenFile('System:AppData/DBTests/redshift/data10,26.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: RedshiftStringTypePatternDoubleOpMoreEq
--connection: RedshiftDBTests
--input: string some_number = ">=975" {pattern: double}
--test: DBTests:expectTable(RedshiftStringTypePatternDoubleOpMoreEq(), OpenFile('System:AppData/DBTests/redshift/data10,26.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: RedshiftStringTypePatternDoubleOpLess
--connection: RedshiftDBTests
--input: string some_number = "<20" {pattern: double}
--test: DBTests:expectTable(RedshiftStringTypePatternDoubleOpLess(), OpenFile('System:AppData/DBTests/redshift/data5.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: RedshiftStringTypePatternDoubleOpLessEq
--connection: RedshiftDBTests
--input: string some_number = "<=20" {pattern: double}
--test: DBTests:expectTable(RedshiftStringTypePatternDoubleOpLessEq(), OpenFile('System:AppData/DBTests/redshift/data5.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: RedshiftStringTypePatternStringOpContains
--connection: RedshiftDBTests
--input: string first_name = "contains Z" {pattern: string}
--test: DBTests:expectTable(RedshiftStringTypePatternStringOpContains(first_name = 'contains Z'), OpenFile('System:AppData/DBTests/redshift/data25.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: RedshiftStringTypePatternStringOpStartsWith
--connection: RedshiftDBTests
--input: string first_name = "starts with W" {pattern: string}
--test: DBTests:expectTable(RedshiftStringTypePatternStringOpStartsWith(first_name='starts with W'), OpenFile('System:AppData/DBTests/redshift/data23.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: RedshiftStringTypePatternStringOpEndsWith
--connection: RedshiftDBTests
--input: string first_name = "ends with y" {pattern: string}
--test: DBTests:expectTable(RedshiftStringTypePatternStringOpEndsWith(first_name = 'ends with y'), OpenFile('System:AppData/DBTests/redshift/data6,23,25.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: RedshiftStringTypePatternStringOpIn
--connection: RedshiftDBTests
--input: string country = "in (Poland, Brazil)" {pattern: string}
--test: DBTests:expectTable(RedshiftStringTypePatternStringOpIn(), OpenFile('System:AppData/DBTests/redshift/data2,5,20.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @country(country);
--end

--name: RedshiftStringTypePatternStringOpRegex
--connection: RedshiftDBTests
--input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
--test: DBTests:expectTable(RedshiftStringTypePatternStringOpRegex(email = 'regex ^([A-Za-z0-9_]+@google.com.au)$'), OpenFile('System:AppData/DBTests/redshift/data9.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @email(email);
--end

--name: RedshiftPatternsAllParams
--connection: RedshiftDBTests
--input: string first_name = "starts with p" {pattern: string}
--input: string id = ">1" {pattern :int}
--input: bool bool = false
--input: string email = "contains com" {pattern: string}
--input: string some_number = ">20" {pattern: double}
--input: string country = "in (Indonesia)" {pattern: string}
--input: string date = "before 1/1/2022" {pattern: datetime}
--test: DBTests:expectTable(RedshiftPatternsAllParams(first_name = "starts with p", id = ">1", false, email = "contains com", some_number = ">20", country = "in (Indonesia)", date = "before 1/1/2022"), OpenFile("System:AppData/DBTests/redshift/data13.d42")) //skip: GROK-12507
SELECT * FROM mock_data
WHERE @first_name(first_name)
  AND @id(id)
           AND bool = @bool
           AND @email(email)
           AND @some_number(some_number)
           AND @country(country)
           AND @date(date);
--end
