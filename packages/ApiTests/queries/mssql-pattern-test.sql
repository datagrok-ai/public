--name: MSSQLPatternsAll
--connection: MSSQLApiTests
--test: ApiTests:expectTable(MSSQLPatternsAll(), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data1-30.d42'))
SELECT * FROM mock_data;
--end

-- INT PATTERN

--name: MSSQLIntTypePatternNone
--connection: MSSQLApiTests
--input: int id = 20
--test: ApiTests:expectTable(MSSQLIntTypePatternNone(), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data20.d42'))
SELECT * FROM mock_data WHERE id = @id;
--end

--name: MSSQLStringTypeIntPatternOpMore
--connection: MSSQLApiTests
--input: string id = ">28" {pattern: int}
--test: ApiTests:expectTable(MSSQLStringTypeIntPatternOpMore(), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: MSSQLStringTypeIntPatternOpMoreEq
--connection: MSSQLApiTests
--input: string id = ">=29" {pattern: int}
--test: ApiTests:expectTable(MSSQLStringTypeIntPatternOpMoreEq(), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: MSSQLStringTypeIntPatternOpLessEq
--connection: MSSQLApiTests
--input: string id = "<=1" {pattern: int}
--test: ApiTests:expectTable(MSSQLStringTypeIntPatternOpLessEq(), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data1.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: MSSQLStringTypeIntPatternOpLess
--connection: MSSQLApiTests
--input: string id = "<2" {pattern: int}
--test: ApiTests:expectTable(MSSQLStringTypeIntPatternOpLess(), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data1.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: MSSQLStringTypeIntPatternOpIn
--connection: MSSQLApiTests
--input: string id = "in(29, 30)" {pattern: int}
--test: ApiTests:expectTable(MSSQLStringTypeIntPatternOpIn(), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: MSSQLStringTypeIntPatternOpNotIn
--connection: MSSQLApiTests
--input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
--test: ApiTests:expectTable(MSSQLStringTypeIntPatternOpNotIn(), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data1-20.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: MSSQLStringTypeIntPatternOpMinMax
--connection: MSSQLApiTests
--input: string id = "min-max 29-30" {pattern: int}
--test: ApiTests:expectTable(MSSQLStringTypeIntPatternOpMinMax(), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: MSSQLStringTypeIntPatternOpNotEq
--connection: MSSQLApiTests
--input: string id = "!=1" {pattern: int}
--test: ApiTests:expectTable(MSSQLStringTypeIntPatternOpNotEq(), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data2-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end


--DOUBLE PATTERN

--name: MSSQLDoubleTypePatternNone
--connection: MSSQLApiTests
--input: double some_number = 510.32
--test: ApiTests:expectTable(MSSQLDoubleTypePatternNone(), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data1.d42'))
SELECT * FROM mock_data WHERE some_number = @some_number;
--end

--name: MSSQLStringTypePatternDoubleOpMore
--connection: MSSQLApiTests
--input: string some_number = ">975" {pattern: double}
--test: ApiTests:expectTable(MSSQLStringTypePatternDoubleOpMore(), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data10,26.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: MSSQLStringTypePatternDoubleOpMoreEq
--connection: MSSQLApiTests
--input: string some_number = ">=975" {pattern: double}
--test: ApiTests:expectTable(MSSQLStringTypePatternDoubleOpMoreEq(), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data10,26.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: MSSQLStringTypePatternDoubleOpLess
--connection: MSSQLApiTests
--input: string some_number = "<20" {pattern: double}
--test: ApiTests:expectTable(MSSQLStringTypePatternDoubleOpLess(), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data5.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: MSSQLStringTypePatternDoubleOpLessEq
--connection: MSSQLApiTests
--input: string some_number = "<=20" {pattern: double}
--test: ApiTests:expectTable(MSSQLStringTypePatternDoubleOpLessEq(), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data5.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--CHOICES - should be used for end-to-end tests

--name: MSSQLByStringChoices
--input: string country = "France" {choices: ["France", "China", "USA", "Finland"]}
SELECT * FROM mock_data WHERE country = @country;
--end

--name: MSSQLByStringChoices
--input: string country = "France" {choices: Query("SELECT DISTINCT country FROM mock_data")}
SELECT * FROM mock_data WHERE country = @country;
--end

--STRING PATTERN

--name: MSSQLStringTypePatternStringOpContains
--connection: MSSQLApiTests
--input: string first_name = "contains Z" {pattern: string}
--test: ApiTests:expectTable(MSSQLStringTypePatternStringOpContains(), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data25.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: MSSQLStringTypePatternStringOpStartsWith
--connection: MSSQLApiTests
--input: string first_name = "starts with W" {pattern: string}
--test: ApiTests:expectTable(MSSQLStringTypePatternStringOpStartsWith(first_name='starts with W'), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data23.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: MSSQLStringTypePatternStringOpEndsWith
--connection: MSSQLApiTests
--input: string first_name = "ends with y" {pattern: string}
--test: ApiTests:expectTable(MSSQLStringTypePatternStringOpEndsWith(first_name = 'ends with y'), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data6,23,25.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: MSSQLStringTypePatternStringOpIn
--connection: MSSQLApiTests
--input: string country = "in (Poland, Brazil)" {pattern: string}
--test: ApiTests:expectTable(MSSQLStringTypePatternStringOpIn(), OpenFile('System:AppData/ApiTests/datasets/tests/mssql/data2,5,20.d42'))
SELECT * FROM mock_data WHERE @country(country);
--end

--name: MSSQLPatternsAllParams
--connection: MSSQLApiTests
--input: string first_name = "starts with p" {pattern: string}
--input: string id = ">1" {pattern :int}
--input: string email = "contains com" {pattern: string}
--input: string some_number = ">20" {pattern: double}
--input: string country = "in (Indonesia)" {pattern: string}
--input: string date = "before 1/1/2022" {pattern: datetime}
--test: ApiTests:expectTable(MSSQLPatternsAllParams(first_name = "starts with p", id = ">1", email = "contains com", some_number = ">20", country = "in (Indonesia)", date = "before 1/1/2022"), OpenFile("System:AppData/ApiTests/datasets/tests/mssql/data13.d42"))
SELECT * FROM mock_data WHERE @first_name(first_name)
AND @id(id) AND @email(email) AND @some_number(some_number) AND @country(country) AND @date(date);
--end
