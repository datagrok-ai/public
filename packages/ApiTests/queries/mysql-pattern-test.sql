-- name: MySqlPatternsAll
-- connection: MySqlApiTests
-- test: ApiTests:expectTable(MySqlPatternsAll(), OpenFile('System:AppData/ApiTests/datasets/test/postgresql/data1-30.d42'))
SELECT * FROM mock_data;
-- end

-- name: MySqlIntTypePatternNone
-- connection: MySqlApiTests
-- input: int id = 20
-- test: ApiTests:expectTable(MySqlIntTypePatternNone(20), OpenFile('System:AppData/ApiTests/datasets/test/postgresql/data20.d42'))
SELECT * FROM mock_data WHERE id = @id;
-- end

-- name: MySqlStringTypeIntPatternOpMore
-- connection: MySqlApiTests
-- input: string id = ">28" {pattern: int}
-- test: ApiTests:expectTable(MySqlStringTypeIntPatternOpMore(), OpenFile('System:AppData/ApiTests/datasets/test/postgresql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MySqlStringTypeIntPatternOpMoreEq
-- connection: MySqlApiTests
-- input: string id = ">=29" {pattern: int}
-- test: ApiTests:expectTable(MySqlStringTypeIntPatternOpMoreEq(), OpenFile('System:AppData/ApiTests/datasets/test/mysql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MySqlStringTypeIntPatternOpLessEq
-- connection: MySqlApiTests
-- input: string id = "<=1" {pattern: int}
-- test: ApiTests:expectTable(MySqlStringTypeIntPatternOpLessEq(), OpenFile('System:AppData/ApiTests/datasets/test/mysql/data1.d42'))
SELECT * FROM mock_data WHERE @id(id)
--
-- end

-- name: MySqlStringTypeIntPatternOpLess
-- connection: MySqlApiTests
-- input: string id = "<2" {pattern: int}
-- test: ApiTests:expectTable(MySqlStringTypeIntPatternOpLess(), OpenFile('System:AppData/ApiTests/datasets/test/mysql/data1.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MySqlStringTypeIntPatternOpIn
-- connection: MySqlApiTests
-- input: string id = "in(29, 30)" {pattern: int}
-- test: ApiTests:expectTable(MySqlStringTypeIntPatternOpIn(), OpenFile('System:AppData/ApiTests/datasets/test/mysql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MySqlStringTypeIntPatternOpNotIn
-- connection: MySqlApiTests
-- input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
-- test: ApiTests:expectTable(MySqlStringTypeIntPatternOpNotIn(), OpenFile('System:AppData/ApiTests/datasets/test/mysql/data1-20.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MySqlStringTypeIntPatternOpMinMax
-- connection: MySqlApiTests
-- input: string id = "min-max 29-30" {pattern: int}
-- test: ApiTests:expectTable(MySqlStringTypeIntPatternOpMinMax(), OpenFile('System:AppData/ApiTests/datasets/test/mysql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MySqlStringTypeIntPatternOpNotEq
-- connection: MySqlApiTests
-- input: string id = "!=1" {pattern: int}
-- test: ApiTests:expectTable(MySqlStringTypeIntPatternOpNotEq(), OpenFile('System:AppData/ApiTests/datasets/test/mysql/data2-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MySqlDoubleTypePatternNone
-- connection: MySqlApiTests
-- input: double some_number = 510.32
-- test: ApiTests:expectTable(MySqlDoubleTypePatternNone(510.32), OpenFile('System:AppData/ApiTests/datasets/test/mysql/data1.d42'))
SELECT * FROM mock_data WHERE some_number = @some_number;
-- end

-- name: MySqlStringTypePatternDoubleOpMore
-- connection: MySqlApiTests
-- input: string some_number = ">975" {pattern: double}
-- test: ApiTests:expectTable(MySqlStringTypePatternDoubleOpMore(), OpenFile('System:AppData/ApiTests/datasets/test/mysql/data10,26.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: MySqlStringTypePatternDoubleOpMoreEq
-- connection: MySqlApiTests
-- input: string some_number = ">=975" {pattern: double}
-- test: ApiTests:expectTable(MySqlStringTypePatternDoubleOpMoreEq(), OpenFile('System:AppData/ApiTests/datasets/test/mysql/data10,26.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: MySqlStringTypePatternDoubleOpLess
-- connection: MySqlApiTests
-- input: string some_number = "<20" {pattern: double}
-- test: ApiTests:expectTable(MySqlStringTypePatternDoubleOpLess(), OpenFile('System:AppData/ApiTests/datasets/test/mysql/data5.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: MySqlStringTypePatternDoubleOpLessEq
-- connection: MySqlApiTests
-- input: string some_number = "<=20" {pattern: double}
-- test: ApiTests:expectTable(MySqlStringTypePatternDoubleOpLessEq(), OpenFile('System:AppData/ApiTests/datasets/test/mysql/data5.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: MySqlStringTypePatternStringOpContains
-- connection: MySqlApiTests
-- input: string first_name = "contains Z" {pattern: string}
-- test: ApiTests:expectTable(MySqlStringTypePatternStringOpContains(first_name = 'contains Z'), OpenFile('System:AppData/ApiTests/datasets/test/mysql/data25.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: MySqlStringTypePatternStringOpStartsWith
-- connection: MySqlApiTests
-- input: string first_name = "starts with W" {pattern: string}
-- test: ApiTests:expectTable(MySqlStringTypePatternStringOpStartsWith(first_name='starts with W'), OpenFile('System:AppData/ApiTests/datasets/test/mysql/data23.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: MySqlStringTypePatternStringOpEndsWith
-- connection: MySqlApiTests
-- input: string first_name = "ends with y" {pattern: string}
-- test: ApiTests:expectTable(MySqlStringTypePatternStringOpEndsWith(first_name = 'ends with y'), OpenFile('System:AppData/ApiTests/datasets/test/mysql/data6,23,25.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: MySqlStringTypePatternStringOpIn
-- connection: MySqlApiTests
-- input: string country = "in (Poland, Brazil)" {pattern: string}
-- test: ApiTests:expectTable(MySqlStringTypePatternStringOpIn(), OpenFile('System:AppData/ApiTests/datasets/test/mysql/data2,5,20.d42'))
SELECT * FROM mock_data WHERE @country(country);
-- end

-- name: MySqlStringTypePatternStringOpRegex
-- connection: MySqlApiTests
-- input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
-- test: ApiTests:expectTable(MySqlStringTypePatternStringOpRegex(email = 'regex ^([A-Za-z0-9_]+@google.com.au)$'), OpenFile('System:AppData/ApiTests/datasets/test/mysql/data9.d42'))
SELECT * FROM mock_data WHERE @email(email);
-- end

-- name: MySqlPatternsAllParams
-- connection: MySqlApiTests
-- input: string first_name = "starts with p" {pattern: string}
-- input: string id = ">1" {pattern :int}
-- input: bool bool = false
-- input: string email = "contains com" {pattern: string}
-- input: string some_number = ">20" {pattern: double}
-- input: string country = "in (Indonesia)" {pattern: string}
-- input: string date = "before 1/1/2022" {pattern: datetime}
-- test: ApiTests:expectTable(MySqlPatternsAllParams(first_name = "starts with p", id = ">1", false, email = "contains com", some_number = ">20", country = "in (Indonesia)", date = "before 1/1/2022"), OpenFile("System:AppData/ApiTests/datasets/test/mysql/data13.d42"))
SELECT * FROM mock_data
WHERE @first_name(first_name)
  AND @id(id)
           AND bool = @bool
           AND @email(email)
           AND @some_number(some_number)
           AND @country(country)
           AND @date(date);
-- end
