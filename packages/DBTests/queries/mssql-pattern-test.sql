--name: MSSQLPatternsAll
--connection: MSSQLDBTests
--test: Dbtests:expectTable(MSSQLPatternsAll(), OpenFile('System:AppData/Dbtests/mssql/data1-30.d42')) // cat: MSSQL
SELECT * FROM mock_data
--end

--name: MSSQLIntTypePatternNone
--connection: MSSQLDBTests
--input: int id = 20
--test: Dbtests:expectTable(MSSQLIntTypePatternNone(20), OpenFile('System:AppData/Dbtests/mssql/data20.d42')) // cat: MSSQL
SELECT * FROM mock_data WHERE id = @id
--end

--name: MSSQLStringTypeIntPatternOpMore
--connection: MSSQLDBTests
--input: string id = ">28" {pattern: int}
--test: Dbtests:expectTable(MSSQLStringTypeIntPatternOpMore('>28'), OpenFile('System:AppData/Dbtests/mssql/data29-30.d42')) // cat: MSSQL
SELECT * FROM mock_data WHERE @id(id)
--end

--name: MSSQLStringTypeIntPatternOpMoreEq
--connection: MSSQLDBTests
--input: string id = ">=29" {pattern: int}
--test: Dbtests:expectTable(MSSQLStringTypeIntPatternOpMoreEq('>=29'), OpenFile('System:AppData/Dbtests/mssql/data29-30.d42')) // cat: MSSQL
SELECT * FROM mock_data WHERE @id(id)
--end

--name: MSSQLStringTypeIntPatternOpLessEq
--connection: MSSQLDBTests
--input: string id = "<=1" {pattern: int}
--test: Dbtests:expectTable(MSSQLStringTypeIntPatternOpLessEq('<=1'), OpenFile('System:AppData/Dbtests/mssql/data1.d42')) // cat: MSSQL
SELECT * FROM mock_data WHERE @id(id)
--end

--name: MSSQLStringTypeIntPatternOpLess
--connection: MSSQLDBTests
--input: string id = "<2" {pattern: int}
--test: Dbtests:expectTable(MSSQLStringTypeIntPatternOpLess('<2'), OpenFile('System:AppData/Dbtests/mssql/data1.d42')) // cat: MSSQL
SELECT * FROM mock_data WHERE @id(id)
--end

--name: MSSQLStringTypeIntPatternOpIn
--connection: MSSQLDBTests
--input: string id = "in(29, 30)" {pattern: int}
--test: Dbtests:expectTable(MSSQLStringTypeIntPatternOpIn(id='in(29, 30)'), OpenFile('System:AppData/Dbtests/mssql/data29-30.d42')) // cat: MSSQL
SELECT * FROM mock_data WHERE @id(id)
--end

--name: MSSQLStringTypeIntPatternOpNotIn
--connection: MSSQLDBTests
--input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
--test: Dbtests:expectTable(MSSQLStringTypeIntPatternOpNotIn(id='not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)'), OpenFile('System:AppData/Dbtests/mssql/data1-20.d42')) // cat: MSSQL
SELECT * FROM mock_data WHERE @id(id)
--end

--name: MSSQLStringTypeIntPatternOpMinMax
--connection: MSSQLDBTests
--input: string id = "min-max 29-30" {pattern: int}
--test: Dbtests:expectTable(MSSQLStringTypeIntPatternOpMinMax(id='min-max 29-30'), OpenFile('System:AppData/Dbtests/mssql/data29-30.d42')) // cat: MSSQL
SELECT * FROM mock_data WHERE @id(id)
--end

--name: MSSQLStringTypeIntPatternOpNotEq
--connection: MSSQLDBTests
--input: string id = "!=1" {pattern: int}
--test: Dbtests:expectTable(MSSQLStringTypeIntPatternOpNotEq('!=1'), OpenFile('System:AppData/Dbtests/mssql/data2-30.d42')) // cat: MSSQL
SELECT * FROM mock_data WHERE @id(id)
--end

--name: MSSQLDoubleTypePatternNone
--connection: MSSQLDBTests
--input: double some_number = 510.32
--test: Dbtests:expectTable(MSSQLDoubleTypePatternNone(510.32), OpenFile('System:AppData/Dbtests/mssql/data1.d42')) // cat: MSSQL
SELECT * FROM mock_data WHERE some_number = @some_number
--end

--name: MSSQLStringTypePatternDoubleOpMore
--connection: MSSQLDBTests
--input: string some_number = ">975" {pattern: double}
--test: Dbtests:expectTable(MSSQLStringTypePatternDoubleOpMore('>975'), OpenFile('System:AppData/Dbtests/mssql/data10,26.d42')) // cat: MSSQL
SELECT * FROM mock_data WHERE @some_number(some_number)
--end

--name: MSSQLStringTypePatternDoubleOpMoreEq
--connection: MSSQLDBTests
--input: string some_number = ">=975" {pattern: double}
--test: Dbtests:expectTable(MSSQLStringTypePatternDoubleOpMoreEq('>=975'), OpenFile('System:AppData/Dbtests/mssql/data10,26.d42')) // cat: MSSQL
SELECT * FROM mock_data WHERE @some_number(some_number)
--end

--name: MSSQLStringTypePatternDoubleOpLess
--connection: MSSQLDBTests
--input: string some_number = "<20" {pattern: double}
--test: Dbtests:expectTable(MSSQLStringTypePatternDoubleOpLess('<20'), OpenFile('System:AppData/Dbtests/mssql/data5.d42')) // cat: MSSQL
SELECT * FROM mock_data WHERE @some_number(some_number)
--end

--name: MSSQLStringTypePatternDoubleOpLessEq
--connection: MSSQLDBTests
--input: string some_number = "<=20" {pattern: double}
--test: Dbtests:expectTable(MSSQLStringTypePatternDoubleOpLessEq('<=20'), OpenFile('System:AppData/Dbtests/mssql/data5.d42')) // cat: MSSQL
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: MSSQLStringTypePatternStringOpContains
--connection: MSSQLDBTests
--input: string first_name = "contains Z" {pattern: string}
--test: Dbtests:expectTable(MSSQLStringTypePatternStringOpContains(first_name = "contains Z"), OpenFile('System:AppData/Dbtests/mssql/data25.d42')) // cat: MSSQL
SELECT * FROM mock_data WHERE @first_name(first_name)
--end

--name: MSSQLStringTypePatternStringOpStartsWith
--connection: MSSQLDBTests
--input: string first_name = "starts with W" {pattern: string}
--test: Dbtests:expectTable(MSSQLStringTypePatternStringOpStartsWith(first_name='starts with W'), OpenFile('System:AppData/Dbtests/mssql/data23.d42')) // cat: MSSQL
SELECT * FROM mock_data WHERE @first_name(first_name)
--end

--name: MSSQLStringTypePatternStringOpEndsWith
--connection: MSSQLDBTests
--input: string first_name = "ends with y" {pattern: string}
--test: Dbtests:expectTable(MSSQLStringTypePatternStringOpEndsWith(first_name = 'ends with y'), OpenFile('System:AppData/Dbtests/mssql/data6,23,25.d42')) // cat: MSSQL
SELECT * FROM mock_data WHERE @first_name(first_name)
--end

--name: MSSQLStringTypePatternStringOpIn
--connection: MSSQLDBTests
--input: string country = "in (Poland, Brazil)" {pattern: string}
--test: Dbtests:expectTable(MSSQLStringTypePatternStringOpIn(country='in (Poland, Brazil)'), OpenFile('System:AppData/Dbtests/mssql/data2,5,20.d42')) // cat: MSSQL
SELECT * FROM mock_data WHERE @country(country)
--end

--name: MSSQLPatternsAllParams
--connection: MSSQLDBTests
--input: string first_name = "starts with p" {pattern: string}
--input: string id = ">1" {pattern :int}
--input: string email = "contains com" {pattern: string}
--input: string some_number = ">20" {pattern: double}
--input: string country = "in (Indonesia)" {pattern: string}
--input: string date = "before 1/1/2022" {pattern: datetime}
--test: Dbtests:expectTable(MSSQLPatternsAllParams(first_name = "starts with p", id = ">1", email = "contains com", some_number = ">20", country = "in (Indonesia)", date = "before 1/1/2022"), OpenFile("System:AppData/Dbtests/mssql/data13.d42")) // cat: MSSQL
SELECT * FROM mock_data WHERE @first_name(first_name)
AND @id(id) AND @email(email) AND @some_number(some_number) AND @country(country) AND @date(date)
--end
