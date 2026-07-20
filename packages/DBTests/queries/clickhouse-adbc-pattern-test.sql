-- name: ClickHouseADBCPatternsAll
-- connection: ClickHouseADBCDBTests
-- test: Dbtests:expectTable(ClickHouseADBCPatternsAll(), OpenFile('System:AppData/Dbtests/common/data1-30.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data;
-- end

-- name: ClickHouseADBCIntTypePatternNone
-- connection: ClickHouseADBCDBTests
-- input: int id = 20
-- test: Dbtests:expectTable(ClickHouseADBCIntTypePatternNone(20), OpenFile('System:AppData/Dbtests/common/data20.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE id = @id;
-- end

-- name: ClickHouseADBCStringTypeIntPatternOpMore
-- connection: ClickHouseADBCDBTests
-- input: string id = ">28" {pattern: int}
-- test: Dbtests:expectTable(ClickHouseADBCStringTypeIntPatternOpMore('>28'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseADBCStringTypeIntPatternOpMoreEq
-- connection: ClickHouseADBCDBTests
-- input: string id = ">=29" {pattern: int}
-- test: Dbtests:expectTable(ClickHouseADBCStringTypeIntPatternOpMoreEq('>=29'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseADBCStringTypeIntPatternOpLessEq
-- connection: ClickHouseADBCDBTests
-- input: string id = "<=1" {pattern: int}
-- test: Dbtests:expectTable(ClickHouseADBCStringTypeIntPatternOpLessEq('<=1'), OpenFile('System:AppData/Dbtests/common/data1.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE @id(id)
--
-- end

-- name: ClickHouseADBCStringTypeIntPatternOpLess
-- connection: ClickHouseADBCDBTests
-- input: string id = "<2" {pattern: int}
-- test: Dbtests:expectTable(ClickHouseADBCStringTypeIntPatternOpLess('<2'), OpenFile('System:AppData/Dbtests/common/data1.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseADBCStringTypeIntPatternOpIn
-- connection: ClickHouseADBCDBTests
-- input: string id = "in(29, 30)" {pattern: int}
-- test: Dbtests:expectTable(ClickHouseADBCStringTypeIntPatternOpIn(id='in(29, 30)'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseADBCStringTypeIntPatternOpNotIn
-- connection: ClickHouseADBCDBTests
-- input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
-- test: Dbtests:expectTable(ClickHouseADBCStringTypeIntPatternOpNotIn(id='not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)'), OpenFile('System:AppData/Dbtests/common/data1-20.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseADBCStringTypeIntPatternOpMinMax
-- connection: ClickHouseADBCDBTests
-- input: string id = "min-max 29-30" {pattern: int}
-- test: Dbtests:expectTable(ClickHouseADBCStringTypeIntPatternOpMinMax(id='min-max 29-30'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseADBCStringTypeIntPatternOpNotEq
-- connection: ClickHouseADBCDBTests
-- input: string id = "!=1" {pattern: int}
-- test: Dbtests:expectTable(ClickHouseADBCStringTypeIntPatternOpNotEq('!=1'), OpenFile('System:AppData/Dbtests/common/data2-30.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseADBCDoubleTypePatternNone
-- connection: ClickHouseADBCDBTests
-- input: double some_number = 510.32
-- test: Dbtests:expectTable(ClickHouseADBCDoubleTypePatternNone(510.32), OpenFile('System:AppData/Dbtests/common/data1.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE some_number = @some_number;
-- end

-- name: ClickHouseADBCStringTypePatternDoubleOpMore
-- connection: ClickHouseADBCDBTests
-- input: string some_number = ">975" {pattern: double}
-- test: Dbtests:expectTable(ClickHouseADBCStringTypePatternDoubleOpMore('>975'), OpenFile('System:AppData/Dbtests/common/data1026.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: ClickHouseADBCStringTypePatternDoubleOpMoreEq
-- connection: ClickHouseADBCDBTests
-- input: string some_number = ">=975" {pattern: double}
-- test: Dbtests:expectTable(ClickHouseADBCStringTypePatternDoubleOpMoreEq('>=975'), OpenFile('System:AppData/Dbtests/common/data1026.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: ClickHouseADBCStringTypePatternDoubleOpLess
-- connection: ClickHouseADBCDBTests
-- input: string some_number = "<20" {pattern: double}
-- test: Dbtests:expectTable(ClickHouseADBCStringTypePatternDoubleOpLess('<20'), OpenFile('System:AppData/Dbtests/common/data5.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: ClickHouseADBCStringTypePatternDoubleOpLessEq
-- connection: ClickHouseADBCDBTests
-- input: string some_number = "<=20" {pattern: double}
-- test: Dbtests:expectTable(ClickHouseADBCStringTypePatternDoubleOpLessEq('<=20'), OpenFile('System:AppData/Dbtests/common/data5.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: ClickHouseADBCStringTypePatternStringOpContains
-- connection: ClickHouseADBCDBTests
-- input: string first_name = "contains Z" {pattern: string}
-- test: Dbtests:expectTable(ClickHouseADBCStringTypePatternStringOpContains(first_name = 'contains Z'), OpenFile('System:AppData/Dbtests/common/data25.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: ClickHouseADBCStringTypePatternStringOpStartsWith
-- connection: ClickHouseADBCDBTests
-- input: string first_name = "starts with W" {pattern: string}
-- test: Dbtests:expectTable(ClickHouseADBCStringTypePatternStringOpStartsWith(first_name='starts with W'), OpenFile('System:AppData/Dbtests/common/data23.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: ClickHouseADBCStringTypePatternStringOpEndsWith
-- connection: ClickHouseADBCDBTests
-- input: string first_name = "ends with y" {pattern: string}
-- test: Dbtests:expectTable(ClickHouseADBCStringTypePatternStringOpEndsWith(first_name = 'ends with y'), OpenFile('System:AppData/Dbtests/common/data62325.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: ClickHouseADBCStringTypePatternStringOpIn
-- connection: ClickHouseADBCDBTests
-- input: string country = "in (Poland, Brazil)" {pattern: string}
-- test: Dbtests:expectTable(ClickHouseADBCStringTypePatternStringOpIn(country='in (Poland, Brazil)'), OpenFile('System:AppData/Dbtests/common/data2520.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE @country(country);
-- end

-- name: ClickHouseADBCStringTypePatternStringOpRegex
-- connection: ClickHouseADBCDBTests
-- input: string email = "regex ^\w+@google.com.au$" {pattern: string}
-- test: Dbtests:expectTable(ClickHouseADBCStringTypePatternStringOpRegex(email = 'regex ^\w+@google.com.au$'), OpenFile('System:AppData/Dbtests/common/data9.d42')) // cat: ClickHouseADBC
SELECT * FROM mock_data WHERE @email(email);
-- end

-- name: ClickHouseADBCPatternsAllParams
-- connection: ClickHouseADBCDBTests
-- input: string first_name = "starts with p" {pattern: string}
-- input: string id = ">1" {pattern :int}
-- input: bool bool = false
-- input: string email = "contains com" {pattern: string}
-- input: string some_number = ">20" {pattern: double}
-- input: string country = "in (Indonesia)" {pattern: string}
-- input: string date = "before 1/1/2022" {pattern: datetime}
-- test: Dbtests:expectTable(ClickHouseADBCPatternsAllParams(first_name = "starts with p", id = ">1", false, email = "contains com", some_number = ">20", country = "in (Indonesia)", date = "before 1/1/2022"), OpenFile("System:AppData/Dbtests/common/data13.d42")) // cat: ClickHouseADBC
SELECT * FROM mock_data
WHERE @first_name(first_name)
  AND @id(id)
           AND bool = @bool
           AND @email(email)
           AND @some_number(some_number)
           AND @country(country)
           AND @date(date);
-- end
