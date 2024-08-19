--name: OraclePatternsAll
--connection: OracleDBTests
--test: Dbtests:expectTable(OraclePatternsAll(), OpenFile('System:AppData/Dbtests/oracle/data1-30.d42')) // cat: Oracle
SELECT * FROM mock_data
--end

--name: OracleIntTypePatternNone
--connection: OracleDBTests
--input: int id = 20
--test: Dbtests:expectTable(OracleIntTypePatternNone(20), OpenFile('System:AppData/Dbtests/oracle/data20.d42')) // cat: Oracle
SELECT * FROM mock_data WHERE id = @id
--end

--name: OracleStringTypeIntPatternOpMore
--connection: OracleDBTests
--input: string id = ">28" {pattern: int}
--test: Dbtests:expectTable(OracleStringTypeIntPatternOpMore('>28'), OpenFile('System:AppData/Dbtests/oracle/data29-30.d42')) // cat: Oracle
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpMoreEq
--connection: OracleDBTests
--input: string id = ">=29" {pattern: int}
--test: Dbtests:expectTable(OracleStringTypeIntPatternOpMoreEq('>=29'), OpenFile('System:AppData/Dbtests/oracle/data29-30.d42')) // cat: Oracle
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpLessEq
--connection: OracleDBTests
--input: string id = "<=1" {pattern: int}
--test: Dbtests:expectTable(OracleStringTypeIntPatternOpLessEq('<=1'), OpenFile('System:AppData/Dbtests/oracle/data1.d42')) // cat: Oracle
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpLess
--connection: OracleDBTests
--input: string id = "<2" {pattern: int}
--test: Dbtests:expectTable(OracleStringTypeIntPatternOpLess('<2'), OpenFile('System:AppData/Dbtests/oracle/data1.d42')) // cat: Oracle
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpIn
--connection: OracleDBTests
--input: string id = "in(29, 30)" {pattern: int}
--test: Dbtests:expectTable(OracleStringTypeIntPatternOpIn(id='in(29, 30)'), OpenFile('System:AppData/Dbtests/oracle/data29-30.d42')) // cat: Oracle
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpNotIn
--connection: OracleDBTests
--input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
--test: Dbtests:expectTable(OracleStringTypeIntPatternOpNotIn(id='not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)'), OpenFile('System:AppData/Dbtests/oracle/data1-20.d42')) // cat: Oracle
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpMinMax
--connection: OracleDBTests
--input: string id = "min-max 29-30" {pattern: int}
--test: Dbtests:expectTable(OracleStringTypeIntPatternOpMinMax(id='min-max 29-30'), OpenFile('System:AppData/Dbtests/oracle/data29-30.d42')) // cat: Oracle
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpNotEq
--connection: OracleDBTests
--input: string id = "!=1" {pattern: int}
--test: Dbtests:expectTable(OracleStringTypeIntPatternOpNotEq('!=1'), OpenFile('System:AppData/Dbtests/oracle/data2-30.d42')) // cat: Oracle
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleDoubleTypePatternNone
--connection: OracleDBTests
--input: double some_number = 510.32
--test: Dbtests:expectTable(OracleDoubleTypePatternNone(510.32), OpenFile('System:AppData/Dbtests/oracle/data1.d42')) // cat: Oracle
SELECT * FROM mock_data WHERE some_number = @some_number
--end

--name: OracleStringTypePatternDoubleOpMore
--connection: OracleDBTests
--input: string some_number = ">975" {pattern: double}
--test: Dbtests:expectTable(OracleStringTypePatternDoubleOpMore('>975'), OpenFile('System:AppData/Dbtests/oracle/data10,26.d42')) // cat: Oracle
SELECT * FROM mock_data WHERE @some_number(some_number)
--end

--name: OracleStringTypePatternDoubleOpMoreEq
--connection: OracleDBTests
--input: string some_number = ">=975" {pattern: double}
--test: Dbtests:expectTable(OracleStringTypePatternDoubleOpMoreEq('>=975'), OpenFile('System:AppData/Dbtests/oracle/data10,26.d42')) // cat: Oracle
SELECT * FROM mock_data WHERE @some_number(some_number)
--end

--name: OracleStringTypePatternDoubleOpLess
--connection: OracleDBTests
--input: string some_number = "<20" {pattern: double}
--test: Dbtests:expectTable(OracleStringTypePatternDoubleOpLess('<20'), OpenFile('System:AppData/Dbtests/oracle/data5.d42')) // cat: Oracle
SELECT * FROM mock_data WHERE @some_number(some_number)
--end

--name: OracleStringTypePatternDoubleOpLessEq
--connection: OracleDBTests
--input: string some_number = "<=20" {pattern: double}
--test: Dbtests:expectTable(OracleStringTypePatternDoubleOpLessEq('<=20'), OpenFile('System:AppData/Dbtests/oracle/data5.d42')) // cat: Oracle
SELECT * FROM mock_data WHERE @some_number(some_number)
--end

--name: OracleStringTypePatternStringOpContains
--connection: OracleDBTests
--input: string first_name = "contains Z" {pattern: string}
--test: Dbtests:expectTable(OracleStringTypePatternStringOpContains(first_name = 'contains Z'), OpenFile('System:AppData/Dbtests/oracle/data25.d42')) // cat: Oracle
SELECT * FROM mock_data WHERE @first_name(first_name)
--end

--name: OracleStringTypePatternStringOpStartsWith
--connection: OracleDBTests
--input: string first_name = "starts with W" {pattern: string}
--test: Dbtests:expectTable(OracleStringTypePatternStringOpStartsWith(first_name='starts with W'), OpenFile('System:AppData/Dbtests/oracle/data23.d42')) // cat: Oracle
SELECT * FROM mock_data WHERE @first_name(first_name)
--end

--name: OracleStringTypePatternStringOpEndsWith
--connection: OracleDBTests
--input: string first_name = "ends with y" {pattern: string}
--test: Dbtests:expectTable(OracleStringTypePatternStringOpEndsWith(first_name = 'ends with y'), OpenFile('System:AppData/Dbtests/oracle/data6,23,25.d42')) // cat: Oracle
SELECT * FROM mock_data WHERE @first_name(first_name)
--end

--name: OracleStringTypePatternStringOpIn
--connection: OracleDBTests
--input: string country = "in (Poland, Brazil)" {pattern: string}
--test: Dbtests:expectTable(OracleStringTypePatternStringOpIn(country='in (Poland, Brazil)'), OpenFile('System:AppData/Dbtests/oracle/data2,5,20.d42')) // cat: Oracle
SELECT * FROM mock_data WHERE @country(country)
--end

--name: OracleStringTypePatternStringOpRegex
--connection: OracleDBTests
--input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
--test: Dbtests:expectTable(OracleStringTypePatternStringOpRegex(email = 'regex ^([A-Za-z0-9_]+@google.com.au)$'), OpenFile('System:AppData/Dbtests/oracle/data9.d42')) // cat: Oracle
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
--test: Dbtests:expectTable(OraclePatternsAllParams(first_name = "starts with p", id = ">1", email = "contains com", some_number = ">20", country = "in (Indonesia)", dat = "before 1/1/2022"), OpenFile("System:AppData/Dbtests/oracle/data13.d42")) // cat: Oracle
SELECT * FROM mock_data WHERE @first_name(first_name)
AND @id(id) AND @email(email) AND @some_number(some_number) AND @country(country) AND @dat(dat)
--end
