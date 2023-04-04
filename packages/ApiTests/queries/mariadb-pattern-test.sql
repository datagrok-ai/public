-- name: MariaDbPatternsAll
-- connection: MariaDbApiTests
-- test: ApiTests:expectTable(MariaDbPatternsAll(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data1-30.d42'))
SELECT * FROM mock_data;
-- end

-- name: MariaDbIntTypePatternNone
-- connection: MariaDbApiTests
-- input: int id = 20
-- test: ApiTests:expectTable(MariaDbIntTypePatternNone(20), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data20.d42'))
SELECT * FROM mock_data WHERE id = @id;
-- end

-- name: MariaDbStringTypeIntPatternOpMore
-- connection: MariaDbApiTests
-- input: string id = ">28" {pattern: int}
-- test: ApiTests:expectTable(MariaDbStringTypeIntPatternOpMore(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MariaDbStringTypeIntPatternOpMoreEq
-- connection: MariaDbApiTests
-- input: string id = ">=29" {pattern: int}
-- test: ApiTests:expectTable(MariaDbStringTypeIntPatternOpMoreEq(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MariaDbStringTypeIntPatternOpLessEq
-- connection: MariaDbApiTests
-- input: string id = "<=1" {pattern: int}
-- test: ApiTests:expectTable(MariaDbStringTypeIntPatternOpLessEq(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data1.d42'))
SELECT * FROM mock_data WHERE @id(id)
--
-- end

-- name: MariaDbStringTypeIntPatternOpLess
-- connection: MariaDbApiTests
-- input: string id = "<2" {pattern: int}
-- test: ApiTests:expectTable(MariaDbStringTypeIntPatternOpLess(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data1.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MariaDbStringTypeIntPatternOpIn
-- connection: MariaDbApiTests
-- input: string id = "in(29, 30)" {pattern: int}
-- test: ApiTests:expectTable(MariaDbStringTypeIntPatternOpIn(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MariaDbStringTypeIntPatternOpNotIn
-- connection: MariaDbApiTests
-- input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
-- test: ApiTests:expectTable(MariaDbStringTypeIntPatternOpNotIn(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data1-20.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MariaDbStringTypeIntPatternOpMinMax
-- connection: MariaDbApiTests
-- input: string id = "min-max 29-30" {pattern: int}
-- test: ApiTests:expectTable(MariaDbStringTypeIntPatternOpMinMax(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MariaDbStringTypeIntPatternOpNotEq
-- connection: MariaDbApiTests
-- input: string id = "!=1" {pattern: int}
-- test: ApiTests:expectTable(MariaDbStringTypeIntPatternOpNotEq(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data2-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MariaDbDoubleTypePatternNone
-- connection: MariaDbApiTests
-- input: double some_number = 510.32
-- test: ApiTests:expectTable(MariaDbDoubleTypePatternNone(510.32), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data1.d42'))
SELECT * FROM mock_data WHERE some_number = @some_number;
-- end

-- name: MariaDbStringTypePatternDoubleOpMore
-- connection: MariaDbApiTests
-- input: string some_number = ">975" {pattern: double}
-- test: ApiTests:expectTable(MariaDbStringTypePatternDoubleOpMore(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data10,26.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: MariaDbStringTypePatternDoubleOpMoreEq
-- connection: MariaDbApiTests
-- input: string some_number = ">=975" {pattern: double}
-- test: ApiTests:expectTable(MariaDbStringTypePatternDoubleOpMoreEq(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data10,26.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: MariaDbStringTypePatternDoubleOpLess
-- connection: MariaDbApiTests
-- input: string some_number = "<20" {pattern: double}
-- test: ApiTests:expectTable(MariaDbStringTypePatternDoubleOpLess(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data5.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: MariaDbStringTypePatternDoubleOpLessEq
-- connection: MariaDbApiTests
-- input: string some_number = "<=20" {pattern: double}
-- test: ApiTests:expectTable(MariaDbStringTypePatternDoubleOpLessEq(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data5.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: MariaDbStringTypePatternStringOpContains
-- connection: MariaDbApiTests
-- input: string first_name = "contains Z" {pattern: string}
-- test: ApiTests:expectTable(MariaDbStringTypePatternStringOpContains(first_name = 'contains Z'), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data25.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: MariaDbStringTypePatternStringOpStartsWith
-- connection: MariaDbApiTests
-- input: string first_name = "starts with W" {pattern: string}
-- test: ApiTests:expectTable(MariaDbStringTypePatternStringOpStartsWith(first_name='starts with W'), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data23.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: MariaDbStringTypePatternStringOpEndsWith
-- connection: MariaDbApiTests
-- input: string first_name = "ends with y" {pattern: string}
-- test: ApiTests:expectTable(MariaDbStringTypePatternStringOpEndsWith(first_name = 'ends with y'), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data6,23,25.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: MariaDbStringTypePatternStringOpIn
-- connection: MariaDbApiTests
-- input: string country = "in (Poland, Brazil)" {pattern: string}
-- test: ApiTests:expectTable(MariaDbStringTypePatternStringOpIn(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data2,5,20.d42'))
SELECT * FROM mock_data WHERE @country(country);
-- end

-- name: MariaDbStringTypePatternStringOpRegex
-- connection: MariaDbApiTests
-- input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
-- test: ApiTests:expectTable(MariaDbStringTypePatternStringOpRegex(email = 'regex ^([A-Za-z0-9_]+@google.com.au)$'), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data9.d42'))
SELECT * FROM mock_data WHERE @email(email);
-- end

-- name: MariaDbPatternsAllParams
-- connection: MariaDbApiTests
-- input: string first_name = "starts with p" {pattern: string}
-- input: string id = ">1" {pattern :int}
-- input: bool bool = false
-- input: string email = "contains com" {pattern: string}
-- input: string some_number = ">20" {pattern: double}
-- input: string country = "in (Indonesia)" {pattern: string}
-- input: string date = "before 1/1/2022" {pattern: datetime}
-- test: ApiTests:expectTable(MariaDbPatternsAllParams(first_name = "starts with p", id = ">1", false, email = "contains com", some_number = ">20", country = "in (Indonesia)", date = "before 1/1/2022"), OpenFile("System:AppData/ApiTests/datasets/tests/postgresql/data13.d42"))
SELECT * FROM mock_data
WHERE @first_name(first_name)
  AND @id(id)
           AND bool = @bool
           AND @email(email)
           AND @some_number(some_number)
           AND @country(country)
           AND @date(date);
-- end
