--name: OraclePatternsAll
--connection: OracleDBTests
--test: DBTests:expectTable(OraclePatternsAll(), OpenFile('System:AppData/DBTests/oracle/data1-30.d42'))
SELECT * FROM mock_data
--end

--name: OracleIntTypePatternNone
--connection: OracleDBTests
--input: int id = 20
--test: DBTests:expectTable(OracleIntTypePatternNone(20), OpenFile('System:AppData/DBTests/oracle/data20.d42'))
SELECT * FROM mock_data WHERE id = @id
--end

--name: OracleStringTypeIntPatternOpMore
--connection: OracleDBTests
--input: string id = ">28" {pattern: int}
--test: DBTests:expectTable(OracleStringTypeIntPatternOpMore(), OpenFile('System:AppData/DBTests/oracle/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpMoreEq
--connection: OracleDBTests
--input: string id = ">=29" {pattern: int}
--test: DBTests:expectTable(OracleStringTypeIntPatternOpMoreEq(), OpenFile('System:AppData/DBTests/oracle/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpLessEq
--connection: OracleDBTests
--input: string id = "<=1" {pattern: int}
--test: DBTests:expectTable(OracleStringTypeIntPatternOpLessEq(), OpenFile('System:AppData/DBTests/oracle/data1.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpLess
--connection: OracleDBTests
--input: string id = "<2" {pattern: int}
--test: DBTests:expectTable(OracleStringTypeIntPatternOpLess(), OpenFile('System:AppData/DBTests/oracle/data1.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpIn
--connection: OracleDBTests
--input: string id = "in(29, 30)" {pattern: int}
--test: DBTests:expectTable(OracleStringTypeIntPatternOpIn(), OpenFile('System:AppData/DBTests/oracle/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpNotIn
--connection: OracleDBTests
--input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
--test: DBTests:expectTable(OracleStringTypeIntPatternOpNotIn(), OpenFile('System:AppData/DBTests/oracle/data1-20.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpMinMax
--connection: OracleDBTests
--input: string id = "min-max 29-30" {pattern: int}
--test: DBTests:expectTable(OracleStringTypeIntPatternOpMinMax(), OpenFile('System:AppData/DBTests/oracle/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpNotEq
--connection: OracleDBTests
--input: string id = "!=1" {pattern: int}
--test: DBTests:expectTable(OracleStringTypeIntPatternOpNotEq(), OpenFile('System:AppData/DBTests/oracle/data2-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleDoubleTypePatternNone
--connection: OracleDBTests
--input: double some_number = 510.32
--test: DBTests:expectTable(OracleDoubleTypePatternNone(510.32), OpenFile('System:AppData/DBTests/oracle/data1.d42'))
SELECT * FROM mock_data WHERE some_number = @some_number
--end

--name: OracleStringTypePatternDoubleOpMore
--connection: OracleDBTests
--input: string some_number = ">975" {pattern: double}
--test: DBTests:expectTable(OracleStringTypePatternDoubleOpMore(), OpenFile('System:AppData/DBTests/oracle/data10,26.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number)
--end

--name: OracleStringTypePatternDoubleOpMoreEq
--connection: OracleDBTests
--input: string some_number = ">=975" {pattern: double}
--test: DBTests:expectTable(OracleStringTypePatternDoubleOpMoreEq(), OpenFile('System:AppData/DBTests/oracle/data10,26.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number)
--end

--name: OracleStringTypePatternDoubleOpLess
--connection: OracleDBTests
--input: string some_number = "<20" {pattern: double}
--test: DBTests:expectTable(OracleStringTypePatternDoubleOpLess(), OpenFile('System:AppData/DBTests/oracle/data5.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number)
--end

--name: OracleStringTypePatternDoubleOpLessEq
--connection: OracleDBTests
--input: string some_number = "<=20" {pattern: double}
--test: DBTests:expectTable(OracleStringTypePatternDoubleOpLessEq(), OpenFile('System:AppData/DBTests/oracle/data5.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number)
--end

--name: OracleStringTypePatternStringOpContains
--connection: OracleDBTests
--input: string first_name = "contains Z" {pattern: string}
--test: DBTests:expectTable(OracleStringTypePatternStringOpContains(first_name = 'contains Z'), OpenFile('System:AppData/DBTests/oracle/data25.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name)
--end

--name: OracleStringTypePatternStringOpStartsWith
--connection: OracleDBTests
--input: string first_name = "starts with W" {pattern: string}
--test: DBTests:expectTable(OracleStringTypePatternStringOpStartsWith(first_name='starts with W'), OpenFile('System:AppData/DBTests/oracle/data23.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name)
--end

--name: OracleStringTypePatternStringOpEndsWith
--connection: OracleDBTests
--input: string first_name = "ends with y" {pattern: string}
--test: DBTests:expectTable(OracleStringTypePatternStringOpEndsWith(first_name = 'ends with y'), OpenFile('System:AppData/DBTests/oracle/data6,23,25.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name)
--end

--name: OracleStringTypePatternStringOpIn
--connection: OracleDBTests
--input: string country = "in (Poland, Brazil)" {pattern: string}
--test: DBTests:expectTable(OracleStringTypePatternStringOpIn(), OpenFile('System:AppData/DBTests/oracle/data2,5,20.d42'))
SELECT * FROM mock_data WHERE @country(country)
--end

--name: OracleStringTypePatternStringOpRegex
--connection: OracleDBTests
--input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
--test: DBTests:expectTable(OracleStringTypePatternStringOpRegex(email = 'regex ^([A-Za-z0-9_]+@google.com.au)$'), OpenFile('System:AppData/DBTests/oracle/data9.d42'))
SELECT * FROM mock_data WHERE @email(email)
--end

--name: OraclePatternsAllParams
--connection: OracleDBTests
--input: string first_name = "starts with p" {pattern: string}
--input: string id = ">1" {pattern :int}
--input: string email = "contains com" {pattern: string}
--input: string some_number = ">20" {pattern: double}
--input: string country = "in (Indonesia)" {pattern: string}
--input: string dat = "before 1/1/2022" {pattern: datetime}
--test: DBTests:expectTable(OraclePatternsAllParams(first_name = "starts with p", id = ">1", email = "contains com", some_number = ">20", country = "in (Indonesia)", dat = "before 1/1/2022"), OpenFile("System:AppData/DBTests/oracle/data13.d42"))
SELECT * FROM mock_data WHERE @first_name(first_name)
AND @id(id) AND @email(email) AND @some_number(some_number) AND @country(country) AND @dat(dat)
--end
