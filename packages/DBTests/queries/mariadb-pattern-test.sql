-- name: MariaDbPatternsAll
-- connection: MariaDbDBTests
-- test: DBTests:expectTable(MariaDbPatternsAll(), OpenFile('System:AppData/DBTests/postgresql/data1-30.d42'))
SELECT * FROM mock_data;
-- end

-- name: MariaDbIntTypePatternNone
-- connection: MariaDbDBTests
-- input: int id = 20
-- test: DBTests:expectTable(MariaDbIntTypePatternNone(20), OpenFile('System:AppData/DBTests/postgresql/data20.d42'))
SELECT * FROM mock_data WHERE id = @id;
-- end

-- name: MariaDbStringTypeIntPatternOpMore
-- connection: MariaDbDBTests
-- input: string id = ">28" {pattern: int}
-- test: DBTests:expectTable(MariaDbStringTypeIntPatternOpMore(), OpenFile('System:AppData/DBTests/postgresql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MariaDbStringTypeIntPatternOpMoreEq
-- connection: MariaDbDBTests
-- input: string id = ">=29" {pattern: int}
-- test: DBTests:expectTable(MariaDbStringTypeIntPatternOpMoreEq(), OpenFile('System:AppData/DBTests/postgresql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MariaDbStringTypeIntPatternOpLessEq
-- connection: MariaDbDBTests
-- input: string id = "<=1" {pattern: int}
-- test: DBTests:expectTable(MariaDbStringTypeIntPatternOpLessEq(), OpenFile('System:AppData/DBTests/postgresql/data1.d42'))
SELECT * FROM mock_data WHERE @id(id)
--
-- end

-- name: MariaDbStringTypeIntPatternOpLess
-- connection: MariaDbDBTests
-- input: string id = "<2" {pattern: int}
-- test: DBTests:expectTable(MariaDbStringTypeIntPatternOpLess(), OpenFile('System:AppData/DBTests/postgresql/data1.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MariaDbStringTypeIntPatternOpIn
-- connection: MariaDbDBTests
-- input: string id = "in(29, 30)" {pattern: int}
-- test: DBTests:expectTable(MariaDbStringTypeIntPatternOpIn(), OpenFile('System:AppData/DBTests/postgresql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MariaDbStringTypeIntPatternOpNotIn
-- connection: MariaDbDBTests
-- input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
-- test: DBTests:expectTable(MariaDbStringTypeIntPatternOpNotIn(), OpenFile('System:AppData/DBTests/postgresql/data1-20.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MariaDbStringTypeIntPatternOpMinMax
-- connection: MariaDbDBTests
-- input: string id = "min-max 29-30" {pattern: int}
-- test: DBTests:expectTable(MariaDbStringTypeIntPatternOpMinMax(), OpenFile('System:AppData/DBTests/postgresql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MariaDbStringTypeIntPatternOpNotEq
-- connection: MariaDbDBTests
-- input: string id = "!=1" {pattern: int}
-- test: DBTests:expectTable(MariaDbStringTypeIntPatternOpNotEq(), OpenFile('System:AppData/DBTests/postgresql/data2-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MariaDbDoubleTypePatternNone
-- connection: MariaDbDBTests
-- input: double some_number = 510.32
-- test: DBTests:expectTable(MariaDbDoubleTypePatternNone(510.32), OpenFile('System:AppData/DBTests/postgresql/data1.d42'))
SELECT * FROM mock_data WHERE some_number = @some_number;
-- end

-- name: MariaDbStringTypePatternDoubleOpMore
-- connection: MariaDbDBTests
-- input: string some_number = ">975" {pattern: double}
-- test: DBTests:expectTable(MariaDbStringTypePatternDoubleOpMore(), OpenFile('System:AppData/DBTests/postgresql/data10,26.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: MariaDbStringTypePatternDoubleOpMoreEq
-- connection: MariaDbDBTests
-- input: string some_number = ">=975" {pattern: double}
-- test: DBTests:expectTable(MariaDbStringTypePatternDoubleOpMoreEq(), OpenFile('System:AppData/DBTests/postgresql/data10,26.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: MariaDbStringTypePatternDoubleOpLess
-- connection: MariaDbDBTests
-- input: string some_number = "<20" {pattern: double}
-- test: DBTests:expectTable(MariaDbStringTypePatternDoubleOpLess(), OpenFile('System:AppData/DBTests/postgresql/data5.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: MariaDbStringTypePatternDoubleOpLessEq
-- connection: MariaDbDBTests
-- input: string some_number = "<=20" {pattern: double}
-- test: DBTests:expectTable(MariaDbStringTypePatternDoubleOpLessEq(), OpenFile('System:AppData/DBTests/postgresql/data5.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: MariaDbStringTypePatternStringOpContains
-- connection: MariaDbDBTests
-- input: string first_name = "contains Z" {pattern: string}
-- test: DBTests:expectTable(MariaDbStringTypePatternStringOpContains(first_name = 'contains Z'), OpenFile('System:AppData/DBTests/postgresql/data25.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: MariaDbStringTypePatternStringOpStartsWith
-- connection: MariaDbDBTests
-- input: string first_name = "starts with W" {pattern: string}
-- test: DBTests:expectTable(MariaDbStringTypePatternStringOpStartsWith(first_name='starts with W'), OpenFile('System:AppData/DBTests/postgresql/data23.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: MariaDbStringTypePatternStringOpEndsWith
-- connection: MariaDbDBTests
-- input: string first_name = "ends with y" {pattern: string}
-- test: DBTests:expectTable(MariaDbStringTypePatternStringOpEndsWith(first_name = 'ends with y'), OpenFile('System:AppData/DBTests/postgresql/data6,23,25.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: MariaDbStringTypePatternStringOpIn
-- connection: MariaDbDBTests
-- input: string country = "in (Poland, Brazil)" {pattern: string}
-- test: DBTests:expectTable(MariaDbStringTypePatternStringOpIn(), OpenFile('System:AppData/DBTests/postgresql/data2,5,20.d42'))
SELECT * FROM mock_data WHERE @country(country);
-- end

-- name: MariaDbStringTypePatternStringOpRegex
-- connection: MariaDbDBTests
-- input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
-- test: DBTests:expectTable(MariaDbStringTypePatternStringOpRegex(email = 'regex ^([A-Za-z0-9_]+@google.com.au)$'), OpenFile('System:AppData/DBTests/postgresql/data9.d42'))
SELECT * FROM mock_data WHERE @email(email);
-- end

-- name: MariaDbPatternsAllParams
-- connection: MariaDbDBTests
-- input: string first_name = "starts with p" {pattern: string}
-- input: string id = ">1" {pattern :int}
-- input: bool bool = false
-- input: string email = "contains com" {pattern: string}
-- input: string some_number = ">20" {pattern: double}
-- input: string country = "in (Indonesia)" {pattern: string}
-- input: string date = "before 1/1/2022" {pattern: datetime}
-- test: DBTests:expectTable(MariaDbPatternsAllParams(first_name = "starts with p", id = ">1", false, email = "contains com", some_number = ">20", country = "in (Indonesia)", date = "before 1/1/2022"), OpenFile("System:AppData/DBTests/postgresql/data13.d42"))
SELECT * FROM mock_data
WHERE @first_name(first_name)
  AND @id(id)
           AND bool = @bool
           AND @email(email)
           AND @some_number(some_number)
           AND @country(country)
           AND @date(date);
-- end
