--name: OraclePatternsAll
--connection: OracleApiTests
--test: ApiTests:expectTable(OraclePatternsAll(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data1-30.d42'))
SELECT * FROM mock_data;
--end

-- INT PATTERN

--name: OracleIntTypePatternNone
--connection: OracleApiTests
--input: int id = 20
--test: ApiTests:expectTable(OracleIntTypePatternNone(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data20.d42'))
SELECT * FROM mock_data WHERE id = @id;
--end

--name: OracleStringTypeIntPatternOpMore
--connection: OracleApiTests
--input: string id = ">28" {pattern: int}
--test: ApiTests:expectTable(OracleStringTypeIntPatternOpMore(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpMoreEq
--connection: OracleApiTests
--input: string id = ">=29" {pattern: int}
--test: ApiTests:expectTable(OracleStringTypeIntPatternOpMoreEq(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpLessEq
--connection: OracleApiTests
--input: string id = "<=1" {pattern: int}
--test: ApiTests:expectTable(OracleStringTypeIntPatternOpLessEq(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data1.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpLess
--connection: OracleApiTests
--input: string id = "<2" {pattern: int}
--test: ApiTests:expectTable(OracleStringTypeIntPatternOpLess(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data1.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpIn
--connection: OracleApiTests
--input: string id = "in(29, 30)" {pattern: int}
--test: ApiTests:expectTable(OracleStringTypeIntPatternOpIn(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpNotIn
--connection: OracleApiTests
--input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
--test: ApiTests:expectTable(OracleStringTypeIntPatternOpNotIn(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data1-20.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpMinMax
--connection: OracleApiTests
--input: string id = "min-max 29-30" {pattern: int}
--test: ApiTests:expectTable(OracleStringTypeIntPatternOpMinMax(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpNotEq
--connection: OracleApiTests
--input: string id = "!=1" {pattern: int}
--test: ApiTests:expectTable(OracleStringTypeIntPatternOpNotEq(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data2-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end


--DOUBLE PATTERN

--name: OracleDoubleTypePatternNone
--connection: OracleApiTests
--input: double some_number = 510.32
--test: ApiTests:expectTable(OracleDoubleTypePatternNone(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data1.d42'))
SELECT * FROM mock_data WHERE some_number = @some_number;
--end

--name: OracleStringTypePatternDoubleOpMore
--connection: OracleApiTests
--input: string some_number = ">975" {pattern: double}
--test: ApiTests:expectTable(OracleStringTypePatternDoubleOpMore(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data10,26.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: OracleStringTypePatternDoubleOpMoreEq
--connection: OracleApiTests
--input: string some_number = ">=975" {pattern: double}
--test: ApiTests:expectTable(OracleStringTypePatternDoubleOpMoreEq(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data10,26.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: OracleStringTypePatternDoubleOpLess
--connection: OracleApiTests
--input: string some_number = "<20" {pattern: double}
--test: ApiTests:expectTable(OracleStringTypePatternDoubleOpLess(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data5.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: OracleStringTypePatternDoubleOpLessEq
--connection: OracleApiTests
--input: string some_number = "<=20" {pattern: double}
--test: ApiTests:expectTable(OracleStringTypePatternDoubleOpLessEq(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data5.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--CHOICES - should be used for end-to-end tests

--name: OracleByStringChoices
--input: string country = "France" {choices: ["France", "China", "USA", "Finland"]}
SELECT * FROM mock_data WHERE country = @country;
--end

--name: OracleByStringChoices
--input: string country = "France" {choices: Query("SELECT DISTINCT country FROM mock_data")}
SELECT * FROM mock_data WHERE country = @country;
--end

--STRING PATTERN

--name: OracleStringTypePatternStringOpContains
--connection: OracleApiTests
--input: string first_name = "contains Z" {pattern: string}
--test: ApiTests:expectTable(OracleStringTypePatternStringOpContains(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data25.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: OracleStringTypePatternStringOpStartsWith
--connection: OracleApiTests
--input: string first_name = "starts with W" {pattern: string}
--test: ApiTests:expectTable(OracleStringTypePatternStringOpStartsWith(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data23.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: OracleStringTypePatternStringOpEndsWith
--connection: OracleApiTests
--input: string first_name = "ends with y" {pattern: string}
--test: ApiTests:expectTable(OracleStringTypePatternStringOpEndsWith(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data6,23,25.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: OracleStringTypePatternStringOpIn
--connection: OracleApiTests
--input: string country = "in (Poland, Brazil)" {pattern: string}
--test: ApiTests:expectTable(OracleStringTypePatternStringOpIn(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data2,5,20.d42'))
SELECT * FROM mock_data WHERE @country(country);
--end

--name: OracleStringTypePatternStringOpRegex
--connection: OracleApiTests
--input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
--test: ApiTests:expectTable(OracleStringTypePatternStringOpRegex(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/data9.d42'))
SELECT * FROM mock_data WHERE @email(email);
--end

--name: OraclePatternsAllParams
--connection: OracleApiTests
--input: string first_name = "starts with p" {pattern: string}
--input: string id = ">1" {pattern :int}
--input: string email = "contains com" {pattern: string}
--input: string some_number = ">20" {pattern: double}
--input: string country = "in (Indonesia)" {pattern: string}
--input: string date = "before 1/1/2022" {pattern: datetime}
--test: ApiTests:expectTable(OraclePatternsAllParams(), OpenFile("System:AppData/ApiTests/datasets/tests/oracle/data13.d42"))
SELECT * FROM mock_data WHERE @first_name(first_name)
AND @id(id) AND @email(email) AND @some_number(some_number) AND @country(country) AND @date(date);
--end
