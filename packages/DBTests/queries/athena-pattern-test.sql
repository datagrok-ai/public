--name: AthenaPatternsAll
--connection: AthenaDBTests
--test: Dbtests:expectTable(AthenaPatternsAll(), OpenFile('System:AppData/Dbtests/common/data1-30.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data
ORDER BY id
--end

--name: AthenaIntTypePatternNone
--connection: AthenaDBTests
--input: int id = 20
--test: Dbtests:expectTable(AthenaIntTypePatternNone(20), OpenFile('System:AppData/Dbtests/common/data20.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE id = @id
ORDER BY id
--end

--name: AthenaStringTypeIntPatternOpMore
--connection: AthenaDBTests
--input: string id = ">28" {pattern: int}
--test: Dbtests:expectTable(AthenaStringTypeIntPatternOpMore('>28'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE @id(id)
ORDER BY id
--end

--name: AthenaStringTypeIntPatternOpMoreEq
--connection: AthenaDBTests
--input: string id = ">=29" {pattern: int}
--test: Dbtests:expectTable(AthenaStringTypeIntPatternOpMoreEq('>=29'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE @id(id)
ORDER BY id
--end

--name: AthenaStringTypeIntPatternOpLessEq
--connection: AthenaDBTests
--input: string id = "<=1" {pattern: int}
--test: Dbtests:expectTable(AthenaStringTypeIntPatternOpLessEq('<=1'), OpenFile('System:AppData/Dbtests/common/data1.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE @id(id)
ORDER BY id
--end

--name: AthenaStringTypeIntPatternOpLess
--connection: AthenaDBTests
--input: string id = "<2" {pattern: int}
--test: Dbtests:expectTable(AthenaStringTypeIntPatternOpLess('<2'), OpenFile('System:AppData/Dbtests/common/data1.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE @id(id)
ORDER BY id
--end

--name: AthenaStringTypeIntPatternOpIn
--connection: AthenaDBTests
--input: string id = "in(29, 30)" {pattern: int}
--test: Dbtests:expectTable(AthenaStringTypeIntPatternOpIn(id='in(29, 30)'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE @id(id)
ORDER BY id
--end

--name: AthenaStringTypeIntPatternOpNotIn
--connection: AthenaDBTests
--input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
--test: Dbtests:expectTable(AthenaStringTypeIntPatternOpNotIn(id='not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)'), OpenFile('System:AppData/Dbtests/common/data1-20.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE @id(id)
ORDER BY id
--end

--name: AthenaStringTypeIntPatternOpMinMax
--connection: AthenaDBTests
--input: string id = "min-max 29-30" {pattern: int}
--test: Dbtests:expectTable(AthenaStringTypeIntPatternOpMinMax(id='min-max 29-30'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE @id(id)
ORDER BY id
--end

--name: AthenaStringTypeIntPatternOpNotEq
--connection: AthenaDBTests
--input: string id = "!=1" {pattern: int}
--test: Dbtests:expectTable(AthenaStringTypeIntPatternOpNotEq('!=1'), OpenFile('System:AppData/Dbtests/common/data2-30.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE @id(id)
ORDER BY id
--end

--name: AthenaDoubleTypePatternNone
--connection: AthenaDBTests
--input: double some_number = 510.32
--test: Dbtests:expectTable(AthenaDoubleTypePatternNone(510.32), OpenFile('System:AppData/Dbtests/common/data1.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE some_number = @some_number
ORDER BY id
--end

--name: AthenaStringTypePatternDoubleOpMore
--connection: AthenaDBTests
--input: string some_number = ">975" {pattern: double}
--test: Dbtests:expectTable(AthenaStringTypePatternDoubleOpMore('>975'), OpenFile('System:AppData/Dbtests/common/data1026.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE @some_number(some_number)
ORDER BY id
--end

--name: AthenaStringTypePatternDoubleOpMoreEq
--connection: AthenaDBTests
--input: string some_number = ">=975" {pattern: double}
--test: Dbtests:expectTable(AthenaStringTypePatternDoubleOpMoreEq('>=975'), OpenFile('System:AppData/Dbtests/common/data1026.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE @some_number(some_number)
ORDER BY id
--end

--name: AthenaStringTypePatternDoubleOpLess
--connection: AthenaDBTests
--input: string some_number = "<20" {pattern: double}
--test: Dbtests:expectTable(AthenaStringTypePatternDoubleOpLess('<20'), OpenFile('System:AppData/Dbtests/common/data5.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE @some_number(some_number)
ORDER BY id
--end

--name: AthenaStringTypePatternDoubleOpLessEq
--connection: AthenaDBTests
--input: string some_number = "<=20" {pattern: double}
--test: Dbtests:expectTable(AthenaStringTypePatternDoubleOpLessEq('<=20'), OpenFile('System:AppData/Dbtests/common/data5.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE @some_number(some_number)
ORDER BY id
--end

--name: AthenaStringTypePatternStringOpContains
--connection: AthenaDBTests
--input: string first_name = "contains Z" {pattern: string}
--test: Dbtests:expectTable(AthenaStringTypePatternStringOpContains(first_name = 'contains Z'), OpenFile('System:AppData/Dbtests/common/data25.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE @first_name(first_name)
ORDER BY id
--end

--name: AthenaStringTypePatternStringOpStartsWith
--connection: AthenaDBTests
--input: string first_name = "starts with W" {pattern: string}
--test: Dbtests:expectTable(AthenaStringTypePatternStringOpStartsWith(first_name='starts with W'), OpenFile('System:AppData/Dbtests/common/data23.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE @first_name(first_name)
ORDER BY id
--end

--name: AthenaStringTypePatternStringOpEndsWith
--connection: AthenaDBTests
--input: string first_name = "ends with y" {pattern: string}
--test: Dbtests:expectTable(AthenaStringTypePatternStringOpEndsWith(first_name = 'ends with y'), OpenFile('System:AppData/Dbtests/common/data62325.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE @first_name(first_name)
ORDER BY id
--end

--name: AthenaStringTypePatternStringOpIn
--connection: AthenaDBTests
--input: string country = "in (Poland, Brazil)" {pattern: string}
--test: Dbtests:expectTable(AthenaStringTypePatternStringOpIn(country='in (Poland, Brazil)'), OpenFile('System:AppData/Dbtests/common/data2520.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE @country(country)
ORDER BY id
--end

--name: AthenaStringTypePatternStringOpRegex
--connection: AthenaDBTests
--input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
--test: Dbtests:expectTable(AthenaStringTypePatternStringOpRegex(email = 'regex ^([A-Za-z0-9_]+@google.com.au)$'), OpenFile('System:AppData/Dbtests/common/data9.d42')) // cat: Athena, timeout: 150000
SELECT * FROM mock_data WHERE @email(email)
ORDER BY id
--end

--name: AthenaPatternsAllParams
--connection: AthenaDBTests
--input: string first_name = "starts with p" {pattern: string}
--input: string id = ">1" {pattern :int}
--input: bool bool = false
--input: string email = "contains com" {pattern: string}
--input: string some_number = ">20" {pattern: double}
--input: string country = "in (Indonesia)" {pattern: string}
--input: string date = "before 1/1/2022" {pattern: datetime}
--test: Dbtests:expectTable(AthenaPatternsAllParams(first_name = "starts with p", id = ">1", false, email = "contains com", some_number = ">20", country = "in (Indonesia)", date = "before 1/1/2022"), OpenFile("System:AppData/Dbtests/common/data13.d42")) // cat: Athena, timeout: 150000
SELECT * FROM mock_data
WHERE @first_name(first_name)
  AND @id(id)
           AND bool = @bool
           AND @email(email)
           AND @some_number(some_number)
           AND @country(country)
           AND @date(date)
ORDER BY id
--end
