-- name: ClickHousePatternsAll
-- connection: ClickHouseApiTests
-- test: ApiTests:expectTable(ClickHousePatternsAll(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data1-30.d42'))
SELECT * FROM mock_data;
-- end

-- name: ClickHouseIntTypePatternNone
-- connection: ClickHouseApiTests
-- input: int id = 20
-- test: ApiTests:expectTable(ClickHouseIntTypePatternNone(20), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data20.d42'))
SELECT * FROM mock_data WHERE id = @id;
-- end

-- name: ClickHouseStringTypeIntPatternOpMore
-- connection: ClickHouseApiTests
-- input: string id = ">28" {pattern: int}
-- test: ApiTests:expectTable(ClickHouseStringTypeIntPatternOpMore(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseStringTypeIntPatternOpMoreEq
-- connection: ClickHouseApiTests
-- input: string id = ">=29" {pattern: int}
-- test: ApiTests:expectTable(ClickHouseStringTypeIntPatternOpMoreEq(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseStringTypeIntPatternOpLessEq
-- connection: ClickHouseApiTests
-- input: string id = "<=1" {pattern: int}
-- test: ApiTests:expectTable(ClickHouseStringTypeIntPatternOpLessEq(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data1.d42'))
SELECT * FROM mock_data WHERE @id(id)
--
-- end

-- name: ClickHouseStringTypeIntPatternOpLess
-- connection: ClickHouseApiTests
-- input: string id = "<2" {pattern: int}
-- test: ApiTests:expectTable(ClickHouseStringTypeIntPatternOpLess(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data1.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseStringTypeIntPatternOpIn
-- connection: ClickHouseApiTests
-- input: string id = "in(29, 30)" {pattern: int}
-- test: ApiTests:expectTable(ClickHouseStringTypeIntPatternOpIn(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseStringTypeIntPatternOpNotIn
-- connection: ClickHouseApiTests
-- input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
-- test: ApiTests:expectTable(ClickHouseStringTypeIntPatternOpNotIn(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data1-20.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseStringTypeIntPatternOpMinMax
-- connection: ClickHouseApiTests
-- input: string id = "min-max 29-30" {pattern: int}
-- test: ApiTests:expectTable(ClickHouseStringTypeIntPatternOpMinMax(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseStringTypeIntPatternOpNotEq
-- connection: ClickHouseApiTests
-- input: string id = "!=1" {pattern: int}
-- test: ApiTests:expectTable(ClickHouseStringTypeIntPatternOpNotEq(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data2-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseDoubleTypePatternNone
-- connection: ClickHouseApiTests
-- input: double some_number = 510.32
-- test: ApiTests:expectTable(ClickHouseDoubleTypePatternNone(510.32), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data1.d42'))
SELECT * FROM mock_data WHERE some_number = @some_number;
-- end

-- name: ClickHouseStringTypePatternDoubleOpMore
-- connection: ClickHouseApiTests
-- input: string some_number = ">975" {pattern: double}
-- test: ApiTests:expectTable(ClickHouseStringTypePatternDoubleOpMore(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data10,26.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: ClickHouseStringTypePatternDoubleOpMoreEq
-- connection: ClickHouseApiTests
-- input: string some_number = ">=975" {pattern: double}
-- test: ApiTests:expectTable(ClickHouseStringTypePatternDoubleOpMoreEq(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data10,26.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: ClickHouseStringTypePatternDoubleOpLess
-- connection: ClickHouseApiTests
-- input: string some_number = "<20" {pattern: double}
-- test: ApiTests:expectTable(ClickHouseStringTypePatternDoubleOpLess(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data5.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: ClickHouseStringTypePatternDoubleOpLessEq
-- connection: ClickHouseApiTests
-- input: string some_number = "<=20" {pattern: double}
-- test: ApiTests:expectTable(ClickHouseStringTypePatternDoubleOpLessEq(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data5.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: ClickHouseStringTypePatternStringOpContains
-- connection: ClickHouseApiTests
-- input: string first_name = "contains Z" {pattern: string}
-- test: ApiTests:expectTable(ClickHouseStringTypePatternStringOpContains(first_name = 'contains Z'), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data25.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: ClickHouseStringTypePatternStringOpStartsWith
-- connection: ClickHouseApiTests
-- input: string first_name = "starts with W" {pattern: string}
-- test: ApiTests:expectTable(ClickHouseStringTypePatternStringOpStartsWith(first_name='starts with W'), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data23.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: ClickHouseStringTypePatternStringOpEndsWith
-- connection: ClickHouseApiTests
-- input: string first_name = "ends with y" {pattern: string}
-- test: ApiTests:expectTable(ClickHouseStringTypePatternStringOpEndsWith(first_name = 'ends with y'), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data6,23,25.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: ClickHouseStringTypePatternStringOpIn
-- connection: ClickHouseApiTests
-- input: string country = "in (Poland, Brazil)" {pattern: string}
-- test: ApiTests:expectTable(ClickHouseStringTypePatternStringOpIn(), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data2,5,20.d42'))
SELECT * FROM mock_data WHERE @country(country);
-- end

-- name: ClickHouseStringTypePatternStringOpRegex
-- connection: ClickHouseApiTests
-- input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
-- test: ApiTests:expectTable(ClickHouseStringTypePatternStringOpRegex(email = 'regex ^([A-Za-z0-9_]+@google.com.au)$'), OpenFile('System:AppData/ApiTests/datasets/tests/postgresql/data9.d42'))
SELECT * FROM mock_data WHERE @email(email);
-- end

-- name: ClickHousePatternsAllParams
-- connection: ClickHouseApiTests
-- input: string first_name = "starts with p" {pattern: string}
-- input: string id = ">1" {pattern :int}
-- input: bool bool = false
-- input: string email = "contains com" {pattern: string}
-- input: string some_number = ">20" {pattern: double}
-- input: string country = "in (Indonesia)" {pattern: string}
-- input: string date = "before 1/1/2022" {pattern: datetime}
-- test: ApiTests:expectTable(ClickHousePatternsAllParams(first_name = "starts with p", id = ">1", false, email = "contains com", some_number = ">20", country = "in (Indonesia)", date = "before 1/1/2022"), OpenFile("System:AppData/ApiTests/datasets/tests/postgresql/data13.d42"))
SELECT * FROM mock_data
WHERE @first_name(first_name)
  AND @id(id)
           AND bool = @bool
           AND @email(email)
           AND @some_number(some_number)
           AND @country(country)
           AND @date(date);
-- end
