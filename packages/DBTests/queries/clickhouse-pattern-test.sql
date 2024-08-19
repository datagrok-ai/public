-- name: ClickHousePatternsAll
-- connection: ClickHouseDBTests
-- test: Dbtests:expectTable(ClickHousePatternsAll(), OpenFile('System:AppData/Dbtests/common/data1-30.d42')) // cat: ClickHouse
SELECT * FROM mock_data;
-- end

-- name: ClickHouseIntTypePatternNone
-- connection: ClickHouseDBTests
-- input: int id = 20
-- test: Dbtests:expectTable(ClickHouseIntTypePatternNone(20), OpenFile('System:AppData/Dbtests/common/data20.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE id = @id;
-- end

-- name: ClickHouseStringTypeIntPatternOpMore
-- connection: ClickHouseDBTests
-- input: string id = ">28" {pattern: int}
-- test: Dbtests:expectTable(ClickHouseStringTypeIntPatternOpMore('>28'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseStringTypeIntPatternOpMoreEq
-- connection: ClickHouseDBTests
-- input: string id = ">=29" {pattern: int}
-- test: Dbtests:expectTable(ClickHouseStringTypeIntPatternOpMoreEq('>=29'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseStringTypeIntPatternOpLessEq
-- connection: ClickHouseDBTests
-- input: string id = "<=1" {pattern: int}
-- test: Dbtests:expectTable(ClickHouseStringTypeIntPatternOpLessEq('<=1'), OpenFile('System:AppData/Dbtests/common/data1.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE @id(id)
--
-- end

-- name: ClickHouseStringTypeIntPatternOpLess
-- connection: ClickHouseDBTests
-- input: string id = "<2" {pattern: int}
-- test: Dbtests:expectTable(ClickHouseStringTypeIntPatternOpLess('<2'), OpenFile('System:AppData/Dbtests/common/data1.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseStringTypeIntPatternOpIn
-- connection: ClickHouseDBTests
-- input: string id = "in(29, 30)" {pattern: int}
-- test: Dbtests:expectTable(ClickHouseStringTypeIntPatternOpIn(id='in(29, 30)'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseStringTypeIntPatternOpNotIn
-- connection: ClickHouseDBTests
-- input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
-- test: Dbtests:expectTable(ClickHouseStringTypeIntPatternOpNotIn(id='not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)'), OpenFile('System:AppData/Dbtests/common/data1-20.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseStringTypeIntPatternOpMinMax
-- connection: ClickHouseDBTests
-- input: string id = "min-max 29-30" {pattern: int}
-- test: Dbtests:expectTable(ClickHouseStringTypeIntPatternOpMinMax(id='min-max 29-30'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseStringTypeIntPatternOpNotEq
-- connection: ClickHouseDBTests
-- input: string id = "!=1" {pattern: int}
-- test: Dbtests:expectTable(ClickHouseStringTypeIntPatternOpNotEq('!=1'), OpenFile('System:AppData/Dbtests/common/data2-30.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE @id(id)
-- end

-- name: ClickHouseDoubleTypePatternNone
-- connection: ClickHouseDBTests
-- input: double some_number = 510.32
-- test: Dbtests:expectTable(ClickHouseDoubleTypePatternNone(510.32), OpenFile('System:AppData/Dbtests/common/data1.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE some_number = @some_number;
-- end

-- name: ClickHouseStringTypePatternDoubleOpMore
-- connection: ClickHouseDBTests
-- input: string some_number = ">975" {pattern: double}
-- test: Dbtests:expectTable(ClickHouseStringTypePatternDoubleOpMore('>975'), OpenFile('System:AppData/Dbtests/common/data10,26.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: ClickHouseStringTypePatternDoubleOpMoreEq
-- connection: ClickHouseDBTests
-- input: string some_number = ">=975" {pattern: double}
-- test: Dbtests:expectTable(ClickHouseStringTypePatternDoubleOpMoreEq('>=975'), OpenFile('System:AppData/Dbtests/common/data10,26.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: ClickHouseStringTypePatternDoubleOpLess
-- connection: ClickHouseDBTests
-- input: string some_number = "<20" {pattern: double}
-- test: Dbtests:expectTable(ClickHouseStringTypePatternDoubleOpLess('<20'), OpenFile('System:AppData/Dbtests/common/data5.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: ClickHouseStringTypePatternDoubleOpLessEq
-- connection: ClickHouseDBTests
-- input: string some_number = "<=20" {pattern: double}
-- test: Dbtests:expectTable(ClickHouseStringTypePatternDoubleOpLessEq('<=20'), OpenFile('System:AppData/Dbtests/common/data5.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE @some_number(some_number);
-- end

-- name: ClickHouseStringTypePatternStringOpContains
-- connection: ClickHouseDBTests
-- input: string first_name = "contains Z" {pattern: string}
-- test: Dbtests:expectTable(ClickHouseStringTypePatternStringOpContains(first_name = 'contains Z'), OpenFile('System:AppData/Dbtests/common/data25.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: ClickHouseStringTypePatternStringOpStartsWith
-- connection: ClickHouseDBTests
-- input: string first_name = "starts with W" {pattern: string}
-- test: Dbtests:expectTable(ClickHouseStringTypePatternStringOpStartsWith(first_name='starts with W'), OpenFile('System:AppData/Dbtests/common/data23.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: ClickHouseStringTypePatternStringOpEndsWith
-- connection: ClickHouseDBTests
-- input: string first_name = "ends with y" {pattern: string}
-- test: Dbtests:expectTable(ClickHouseStringTypePatternStringOpEndsWith(first_name = 'ends with y'), OpenFile('System:AppData/Dbtests/common/data6,23,25.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE @first_name(first_name);
-- end

-- name: ClickHouseStringTypePatternStringOpIn
-- connection: ClickHouseDBTests
-- input: string country = "in (Poland, Brazil)" {pattern: string}
-- test: Dbtests:expectTable(ClickHouseStringTypePatternStringOpIn(country='in (Poland, Brazil)'), OpenFile('System:AppData/Dbtests/common/data2,5,20.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE @country(country);
-- end

-- name: ClickHouseStringTypePatternStringOpRegex
-- connection: ClickHouseDBTests
-- input: string email = "regex ^\w+@google.com.au$" {pattern: string}
-- test: Dbtests:expectTable(ClickHouseStringTypePatternStringOpRegex(email = 'regex ^\w+@google.com.au$'), OpenFile('System:AppData/Dbtests/common/data9.d42')) // cat: ClickHouse
SELECT * FROM mock_data WHERE @email(email);
-- end

-- name: ClickHousePatternsAllParams
-- connection: ClickHouseDBTests
-- input: string first_name = "starts with p" {pattern: string}
-- input: string id = ">1" {pattern :int}
-- input: bool bool = false
-- input: string email = "contains com" {pattern: string}
-- input: string some_number = ">20" {pattern: double}
-- input: string country = "in (Indonesia)" {pattern: string}
-- input: string date = "before 1/1/2022" {pattern: datetime}
-- test: Dbtests:expectTable(ClickHousePatternsAllParams(first_name = "starts with p", id = ">1", false, email = "contains com", some_number = ">20", country = "in (Indonesia)", date = "before 1/1/2022"), OpenFile("System:AppData/Dbtests/common/data13.d42")) // cat: ClickHouse
SELECT * FROM mock_data
WHERE @first_name(first_name)
  AND @id(id)
           AND bool = @bool
           AND @email(email)
           AND @some_number(some_number)
           AND @country(country)
           AND @date(date);
-- end
