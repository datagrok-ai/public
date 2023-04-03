--name: RedshiftPatternsAll
--connection: RedshiftApiTests
--test: ApiTests:expectTable(RedshiftPatternsAll(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data1-30.d42')) //skip: GROK-12507
SELECT * FROM mock_data
--end

--name: RedshiftIntTypePatternNone
--connection: RedshiftApiTests
--input: int id = 20
--test: ApiTests:expectTable(RedshiftIntTypePatternNone(20), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data20.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE id = @id;
--end

--name: RedshiftStringTypeIntPatternOpMore
--connection: RedshiftApiTests
--input: string id = ">28" {pattern: int}
--test: ApiTests:expectTable(RedshiftStringTypeIntPatternOpMore(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data29-30.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @id(id)
--end

--name: RedshiftStringTypeIntPatternOpMoreEq
--connection: RedshiftApiTests
--input: string id = ">=29" {pattern: int}
--test: ApiTests:expectTable(RedshiftStringTypeIntPatternOpMoreEq(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data29-30.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @id(id)
--end

--name: RedshiftStringTypeIntPatternOpLessEq
--connection: RedshiftApiTests
--input: string id = "<=1" {pattern: int}
--test: ApiTests:expectTable(RedshiftStringTypeIntPatternOpLessEq(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data1.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @id(id)
--end

--name: RedshiftStringTypeIntPatternOpLess
--connection: RedshiftApiTests
--input: string id = "<2" {pattern: int}
--test: ApiTests:expectTable(RedshiftStringTypeIntPatternOpLess(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data1.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @id(id)
--end

--name: RedshiftStringTypeIntPatternOpIn
--connection: RedshiftApiTests
--input: string id = "in(29, 30)" {pattern: int}
--test: ApiTests:expectTable(RedshiftStringTypeIntPatternOpIn(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data29-30.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @id(id)
--end

--name: RedshiftStringTypeIntPatternOpNotIn
--connection: RedshiftApiTests
--input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
--test: ApiTests:expectTable(RedshiftStringTypeIntPatternOpNotIn(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data1-20.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @id(id)
--end

--name: RedshiftStringTypeIntPatternOpMinMax
--connection: RedshiftApiTests
--input: string id = "min-max 29-30" {pattern: int}
--test: ApiTests:expectTable(RedshiftStringTypeIntPatternOpMinMax(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data29-30.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @id(id)
--end

--name: RedshiftStringTypeIntPatternOpNotEq
--connection: RedshiftApiTests
--input: string id = "!=1" {pattern: int}
--test: ApiTests:expectTable(RedshiftStringTypeIntPatternOpNotEq(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data2-30.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @id(id)
--end

--name: RedshiftDoubleTypePatternNone
--connection: RedshiftApiTests
--input: double some_number = 510.32
--test: ApiTests:expectTable(RedshiftDoubleTypePatternNone(510.32), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data1.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE some_number = @some_number;
--end

--name: RedshiftStringTypePatternDoubleOpMore
--connection: RedshiftApiTests
--input: string some_number = ">975" {pattern: double}
--test: ApiTests:expectTable(RedshiftStringTypePatternDoubleOpMore(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data10,26.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: RedshiftStringTypePatternDoubleOpMoreEq
--connection: RedshiftApiTests
--input: string some_number = ">=975" {pattern: double}
--test: ApiTests:expectTable(RedshiftStringTypePatternDoubleOpMoreEq(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data10,26.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: RedshiftStringTypePatternDoubleOpLess
--connection: RedshiftApiTests
--input: string some_number = "<20" {pattern: double}
--test: ApiTests:expectTable(RedshiftStringTypePatternDoubleOpLess(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data5.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: RedshiftStringTypePatternDoubleOpLessEq
--connection: RedshiftApiTests
--input: string some_number = "<=20" {pattern: double}
--test: ApiTests:expectTable(RedshiftStringTypePatternDoubleOpLessEq(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data5.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: RedshiftStringTypePatternStringOpContains
--connection: RedshiftApiTests
--input: string first_name = "contains Z" {pattern: string}
--test: ApiTests:expectTable(RedshiftStringTypePatternStringOpContains(first_name = 'contains Z'), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data25.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: RedshiftStringTypePatternStringOpStartsWith
--connection: RedshiftApiTests
--input: string first_name = "starts with W" {pattern: string}
--test: ApiTests:expectTable(RedshiftStringTypePatternStringOpStartsWith(first_name='starts with W'), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data23.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: RedshiftStringTypePatternStringOpEndsWith
--connection: RedshiftApiTests
--input: string first_name = "ends with y" {pattern: string}
--test: ApiTests:expectTable(RedshiftStringTypePatternStringOpEndsWith(first_name = 'ends with y'), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data6,23,25.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: RedshiftStringTypePatternStringOpIn
--connection: RedshiftApiTests
--input: string country = "in (Poland, Brazil)" {pattern: string}
--test: ApiTests:expectTable(RedshiftStringTypePatternStringOpIn(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data2,5,20.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @country(country);
--end

--name: RedshiftStringTypePatternStringOpRegex
--connection: RedshiftApiTests
--input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
--test: ApiTests:expectTable(RedshiftStringTypePatternStringOpRegex(email = 'regex ^([A-Za-z0-9_]+@google.com.au)$'), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/data9.d42')) //skip: GROK-12507
SELECT * FROM mock_data WHERE @email(email);
--end

--name: RedshiftPatternsAllParams
--connection: RedshiftApiTests
--input: string first_name = "starts with p" {pattern: string}
--input: string id = ">1" {pattern :int}
--input: bool bool = false
--input: string email = "contains com" {pattern: string}
--input: string some_number = ">20" {pattern: double}
--input: string country = "in (Indonesia)" {pattern: string}
--input: string date = "before 1/1/2022" {pattern: datetime}
--test: ApiTests:expectTable(RedshiftPatternsAllParams(first_name = "starts with p", id = ">1", false, email = "contains com", some_number = ">20", country = "in (Indonesia)", date = "before 1/1/2022"), OpenFile("System:AppData/ApiTests/datasets/tests/redshift/data13.d42")) //skip: GROK-12507
SELECT * FROM mock_data
WHERE @first_name(first_name)
  AND @id(id)
           AND bool = @bool
           AND @email(email)
           AND @some_number(some_number)
           AND @country(country)
           AND @date(date);
--end
