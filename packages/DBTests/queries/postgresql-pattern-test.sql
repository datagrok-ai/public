--name: PostgresqlPatternsAll
--connection: PostgreSQLDBTests
--test: Dbtests:expectTable(PostgresqlPatternsAll(), OpenFile('System:AppData/Dbtests/common/data1-30.d42')) // cat: Postgresql
SELECT * FROM mock_data
--end

--name: PostgresqlIntTypePatternNone
--connection: PostgreSQLDBTests
--input: int id = 20
--test: Dbtests:expectTable(PostgresqlIntTypePatternNone(20), OpenFile('System:AppData/Dbtests/common/data20.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE id = @id;
--end

--name: PostgresqlStringTypeIntPatternOpMore
--connection: PostgreSQLDBTests
--input: string id = ">28" {pattern: int}
--test: Dbtests:expectTable(PostgresqlStringTypeIntPatternOpMore('>28'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlStringTypeIntPatternOpMoreEq
--connection: PostgreSQLDBTests
--input: string id = ">=29" {pattern: int}
--test: Dbtests:expectTable(PostgresqlStringTypeIntPatternOpMoreEq('>=29'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlStringTypeIntPatternOpLessEq
--connection: PostgreSQLDBTests
--input: string id = "<=1" {pattern: int}
--test: Dbtests:expectTable(PostgresqlStringTypeIntPatternOpLessEq('<=1'), OpenFile('System:AppData/Dbtests/common/data1.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlStringTypeIntPatternOpLess
--connection: PostgreSQLDBTests
--input: string id = "<2" {pattern: int}
--test: Dbtests:expectTable(PostgresqlStringTypeIntPatternOpLess('<2'), OpenFile('System:AppData/Dbtests/common/data1.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlStringTypeIntPatternOpIn
--connection: PostgreSQLDBTests
--input: string id = "in(29, 30)" {pattern: int}
--test: Dbtests:expectTable(PostgresqlStringTypeIntPatternOpIn(id='in(29, 30)'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlStringTypeIntPatternOpNotIn
--connection: PostgreSQLDBTests
--input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
--test: Dbtests:expectTable(PostgresqlStringTypeIntPatternOpNotIn(id='not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)'), OpenFile('System:AppData/Dbtests/common/data1-20.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlStringTypeIntPatternOpMinMax
--connection: PostgreSQLDBTests
--input: string id = "min-max 29-30" {pattern: int}
--test: Dbtests:expectTable(PostgresqlStringTypeIntPatternOpMinMax(id='min-max 29-30'), OpenFile('System:AppData/Dbtests/common/data29-30.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlStringTypeIntPatternOpNotEq
--connection: PostgreSQLDBTests
--input: string id = "!=1" {pattern: int}
--test: Dbtests:expectTable(PostgresqlStringTypeIntPatternOpNotEq('!=1'), OpenFile('System:AppData/Dbtests/common/data2-30.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlDoubleTypePatternNone
--connection: PostgreSQLDBTests
--input: double some_number = 510.32
--test: Dbtests:expectTable(PostgresqlDoubleTypePatternNone(510.32), OpenFile('System:AppData/Dbtests/common/data1.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE some_number = @some_number;
--end

--name: PostgresqlStringTypePatternDoubleOpMore
--connection: PostgreSQLDBTests
--input: string some_number = ">975" {pattern: double}
--test: Dbtests:expectTable(PostgresqlStringTypePatternDoubleOpMore('>975'), OpenFile('System:AppData/Dbtests/common/data10,26.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: PostgresqlStringTypePatternDoubleOpMoreEq
--connection: PostgreSQLDBTests
--input: string some_number = ">=975" {pattern: double}
--test: Dbtests:expectTable(PostgresqlStringTypePatternDoubleOpMoreEq('>=975'), OpenFile('System:AppData/Dbtests/common/data10,26.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: PostgresqlStringTypePatternDoubleOpLess
--connection: PostgreSQLDBTests
--input: string some_number = "<20" {pattern: double}
--test: Dbtests:expectTable(PostgresqlStringTypePatternDoubleOpLess('<20'), OpenFile('System:AppData/Dbtests/common/data5.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: PostgresqlStringTypePatternDoubleOpLessEq
--connection: PostgreSQLDBTests
--input: string some_number = "<=20" {pattern: double}
--test: Dbtests:expectTable(PostgresqlStringTypePatternDoubleOpLessEq('<=20'), OpenFile('System:AppData/Dbtests/common/data5.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: PostgresqlStringTypePatternStringOpContains
--connection: PostgreSQLDBTests
--input: string first_name = "contains Z" {pattern: string}
--test: Dbtests:expectTable(PostgresqlStringTypePatternStringOpContains(first_name = 'contains Z'), OpenFile('System:AppData/Dbtests/common/data25.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: PostgresqlStringTypePatternStringOpStartsWith
--connection: PostgreSQLDBTests
--input: string first_name = "starts with W" {pattern: string}
--test: Dbtests:expectTable(PostgresqlStringTypePatternStringOpStartsWith(first_name='starts with W'), OpenFile('System:AppData/Dbtests/common/data23.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: PostgresqlStringTypePatternStringOpEndsWith
--connection: PostgreSQLDBTests
--input: string first_name = "ends with y" {pattern: string}
--test: Dbtests:expectTable(PostgresqlStringTypePatternStringOpEndsWith(first_name = 'ends with y'), OpenFile('System:AppData/Dbtests/common/data6,23,25.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: PostgresqlStringTypePatternStringOpIn
--connection: PostgreSQLDBTests
--input: string country = "in (Poland, Brazil)" {pattern: string}
--test: Dbtests:expectTable(PostgresqlStringTypePatternStringOpIn(country='in (Poland, Brazil)'), OpenFile('System:AppData/Dbtests/common/data2,5,20.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE @country(country);
--end

--name: PostgresqlStringTypePatternStringOpRegex
--connection: PostgreSQLDBTests
--input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
--test: Dbtests:expectTable(PostgresqlStringTypePatternStringOpRegex(email = 'regex ^([A-Za-z0-9_]+@google.com.au)$'), OpenFile('System:AppData/Dbtests/common/data9.d42')) // cat: Postgresql
SELECT * FROM mock_data WHERE @email(email);
--end

--name: PostgresqlPatternsAllParams
--connection: PostgreSQLDBTests
--input: string first_name = "starts with p" {pattern: string}
--input: string id = ">1" {pattern :int}
--input: bool bool = false
--input: string email = "contains com" {pattern: string}
--input: string some_number = ">20" {pattern: double}
--input: string country = "in (Indonesia)" {pattern: string}
--input: string date = "before 1/1/2022" {pattern: datetime}
--test: Dbtests:expectTable(PostgresqlPatternsAllParams(first_name = "starts with p", id = ">1", false, email = "contains com", some_number = ">20", country = "in (Indonesia)", date = "before 1/1/2022"), OpenFile("System:AppData/Dbtests/common/data13.d42")) // cat: Postgresql
SELECT * FROM mock_data
WHERE @first_name(first_name)
  AND @id(id)
           AND bool = @bool
           AND @email(email)
           AND @some_number(some_number)
           AND @country(country)
           AND @date(date);
--end
