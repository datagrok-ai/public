-- name: MySqlPatternsAll
-- connection: MySqlDBTests
-- test: Dbtests:expectTable(MySqlPatternsAll(), OpenFile('System:AppData/Dbtests/common/data1-30.d42')) // cat: MySql
SELECT * FROM mock_data;
-- end

-- name: MySqlIntTypePatternNone
-- connection: MySqlDBTests
-- input: int id = 20
-- test: Dbtests:expectTable(MySqlIntTypePatternNone(20), OpenFile('System:AppData/Dbtests/common/data20.d42')) // cat: MySql
SELECT * FROM mock_data WHERE id = @id;
-- end

-- name: MySqlStringTypeIntPatternOpMore
-- connection: MySqlDBTests
-- input: string id = ">28" {pattern: int}
-- test: Dbtests:expectTable(MySqlStringTypeIntPatternOpMore('>28'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: MySql
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MySqlStringTypeIntPatternOpMoreEq
-- connection: MySqlDBTests
-- input: string id = ">=29" {pattern: int}
-- test: Dbtests:expectTable(MySqlStringTypeIntPatternOpMoreEq('>=29'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: MySql
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MySqlStringTypeIntPatternOpLessEq
-- connection: MySqlDBTests
-- input: string id = "<=1" {pattern: int}
-- test: Dbtests:expectTable(MySqlStringTypeIntPatternOpLessEq('<=1'), OpenFile('System:AppData/Dbtests/common/data1.d42')) // cat: MySql
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MySqlStringTypeIntPatternOpLess
-- connection: MySqlDBTests
-- input: string id = "<2" {pattern: int}
-- test: Dbtests:expectTable(MySqlStringTypeIntPatternOpLess('<2'), OpenFile('System:AppData/Dbtests/common/data1.d42')) // cat: MySql
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MySqlStringTypeIntPatternOpIn
-- connection: MySqlDBTests
-- input: string id = "in(29, 30)" {pattern: int}
-- test: Dbtests:expectTable(MySqlStringTypeIntPatternOpIn(id='in(29, 30)'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: MySql
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MySqlStringTypeIntPatternOpNotIn
-- connection: MySqlDBTests
-- input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
-- test: Dbtests:expectTable(MySqlStringTypeIntPatternOpNotIn(id='not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)'), OpenFile('System:AppData/Dbtests/common/data1-20.d42')) // cat: MySql
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MySqlStringTypeIntPatternOpMinMax
-- connection: MySqlDBTests
-- input: string id = "min-max 29-30" {pattern: int}
-- test: Dbtests:expectTable(MySqlStringTypeIntPatternOpMinMax(id='min-max 29-30'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: MySql
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MySqlStringTypeIntPatternOpNotEq
-- connection: MySqlDBTests
-- input: string id = "!=1" {pattern: int}
-- test: Dbtests:expectTable(MySqlStringTypeIntPatternOpNotEq('!=1'), OpenFile('System:AppData/Dbtests/common/data2-30.d42')) // cat: MySql
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: MySqlDoubleTypePatternNone
-- connection: MySqlDBTests
-- input: double some_number = 510.32
-- test: Dbtests:expectTable(MySqlDoubleTypePatternNone(510.32), OpenFile('System:AppData/Dbtests/common/data1.d42')) // cat: MySql
SELECT * FROM mock_data WHERE some_number = @some_number;
-- end

-- name: MySqlStringTypePatternDoubleOpMore
-- connection: MySqlDBTests
-- input: string some_number = ">975" {pattern: double}
-- test: Dbtests:expectTable(MySqlStringTypePatternDoubleOpMore('>975'), OpenFile('System:AppData/Dbtests/common/data10,26.d42')) // cat: MySql
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: MySqlStringTypePatternDoubleOpMoreEq
-- connection: MySqlDBTests
-- input: string some_number = ">=975" {pattern: double}
-- test: Dbtests:expectTable(MySqlStringTypePatternDoubleOpMoreEq('>=975'), OpenFile('System:AppData/Dbtests/common/data10,26.d42')) // cat: MySql
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: MySqlStringTypePatternDoubleOpLess
-- connection: MySqlDBTests
-- input: string some_number = "<20" {pattern: double}
-- test: Dbtests:expectTable(MySqlStringTypePatternDoubleOpLess('<20'), OpenFile('System:AppData/Dbtests/common/data5.d42')) // cat: MySql
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: MySqlStringTypePatternDoubleOpLessEq
-- connection: MySqlDBTests
-- input: string some_number = "<=20" {pattern: double}
-- test: Dbtests:expectTable(MySqlStringTypePatternDoubleOpLessEq('<=20'), OpenFile('System:AppData/Dbtests/common/data5.d42')) // cat: MySql
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: MySqlStringTypePatternStringOpContains
-- connection: MySqlDBTests
-- input: string first_name = "contains Z" {pattern: string}
-- test: Dbtests:expectTable(MySqlStringTypePatternStringOpContains(first_name = 'contains Z'), OpenFile('System:AppData/Dbtests/common/data25.d42')) // cat: MySql
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: MySqlStringTypePatternStringOpStartsWith
-- connection: MySqlDBTests
-- input: string first_name = "starts with W" {pattern: string}
-- test: Dbtests:expectTable(MySqlStringTypePatternStringOpStartsWith(first_name='starts with W'), OpenFile('System:AppData/Dbtests/common/data23.d42')) // cat: MySql
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: MySqlStringTypePatternStringOpEndsWith
-- connection: MySqlDBTests
-- input: string first_name = "ends with y" {pattern: string}
-- test: Dbtests:expectTable(MySqlStringTypePatternStringOpEndsWith(first_name = 'ends with y'), OpenFile('System:AppData/Dbtests/common/data6,23,25.d42')) // cat: MySql
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: MySqlStringTypePatternStringOpIn
-- connection: MySqlDBTests
-- input: string country = "in (Poland, Brazil)" {pattern: string}
-- test: Dbtests:expectTable(MySqlStringTypePatternStringOpIn(country='in (Poland, Brazil)'), OpenFile('System:AppData/Dbtests/common/data2,5,20.d42')) // cat: MySql
SELECT * FROM mock_data WHERE @country(country);
-- end

-- name: MySqlStringTypePatternStringOpRegex
-- connection: MySqlDBTests
-- input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
-- test: Dbtests:expectTable(MySqlStringTypePatternStringOpRegex(email = 'regex ^([A-Za-z0-9_]+@google.com.au)$'), OpenFile('System:AppData/Dbtests/common/data9.d42')) // cat: MySql
SELECT * FROM mock_data WHERE @email(email);
-- end

-- name: MySqlPatternsAllParams
-- connection: MySqlDBTests
-- input: string first_name = "starts with p" {pattern: string}
-- input: string id = ">1" {pattern :int}
-- input: bool bool = false
-- input: string email = "contains com" {pattern: string}
-- input: string some_number = ">20" {pattern: double}
-- input: string country = "in (Indonesia)" {pattern: string}
-- input: string date = "before 1/1/2022" {pattern: datetime}
-- test: Dbtests:expectTable(MySqlPatternsAllParams(first_name = "starts with p", id = ">1", false, email = "contains com", some_number = ">20", country = "in (Indonesia)", date = "before 1/1/2022"), OpenFile("System:AppData/Dbtests/common/data13.d42")) // cat: MySql
SELECT * FROM mock_data
WHERE @first_name(first_name)
  AND @id(id)
           AND bool = @bool
           AND @email(email)
           AND @some_number(some_number)
           AND @country(country)
           AND @date(date);
-- end
