--name: AthenaPatternsAll
--connection: AthenaApiTests
--test: ApiTests:expectTable(AthenaPatternsAll(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data1-30.d42'))
SELECT * FROM mock_data
--end

-- INT PATTERN

--name: AthenaIntTypePatternNone
--connection: AthenaApiTests
--input: int id = 20
--test: ApiTests:expectTable(AthenaIntTypePatternNone(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data20.d42'))
SELECT * FROM mock_data WHERE id = @id;
--end

--name: AthenaStringTypeIntPatternOpMore
--connection: AthenaApiTests
--input: string id = ">28" {pattern: int}
--test: ApiTests:expectTable(AthenaStringTypeIntPatternOpMore(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: AthenaStringTypeIntPatternOpMoreEq
--connection: AthenaApiTests
--input: string id = ">=29" {pattern: int}
--test: ApiTests:expectTable(AthenaStringTypeIntPatternOpMoreEq(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: AthenaStringTypeIntPatternOpLessEq
--connection: AthenaApiTests
--input: string id = "<=1" {pattern: int}
--test: ApiTests:expectTable(AthenaStringTypeIntPatternOpLessEq(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data1.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: AthenaStringTypeIntPatternOpLess
--connection: AthenaApiTests
--input: string id = "<2" {pattern: int}
--test: ApiTests:expectTable(AthenaStringTypeIntPatternOpLess(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data1.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: AthenaStringTypeIntPatternOpIn
--connection: AthenaApiTests
--input: string id = "in(29, 30)" {pattern: int}
--test: ApiTests:expectTable(AthenaStringTypeIntPatternOpIn(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: AthenaStringTypeIntPatternOpNotIn
--connection: AthenaApiTests
--input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
--test: ApiTests:expectTable(AthenaStringTypeIntPatternOpNotIn(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data1-20.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: AthenaStringTypeIntPatternOpMinMax
--connection: AthenaApiTests
--input: string id = "min-max 29-30" {pattern: int}
--test: ApiTests:expectTable(AthenaStringTypeIntPatternOpMinMax(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data29-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: AthenaStringTypeIntPatternOpNotEq
--connection: AthenaApiTests
--input: string id = "!=1" {pattern: int}
--test: ApiTests:expectTable(AthenaStringTypeIntPatternOpNotEq(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data2-30.d42'))
SELECT * FROM mock_data WHERE @id(id)
--end


--DOUBLE PATTERN

--name: AthenaDoubleTypePatternNone
--connection: AthenaApiTests
--input: double some_number = 510.32
--test: ApiTests:expectTable(AthenaDoubleTypePatternNone(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data1.d42'))
SELECT * FROM mock_data WHERE some_number = @some_number;
--end

--name: AthenaStringTypePatternDoubleOpMore
--connection: AthenaApiTests
--input: string some_number = ">975" {pattern: double}
--test: ApiTests:expectTable(AthenaStringTypePatternDoubleOpMore(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data10,26.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: AthenaStringTypePatternDoubleOpMoreEq
--connection: AthenaApiTests
--input: string some_number = ">=975" {pattern: double}
--test: ApiTests:expectTable(AthenaStringTypePatternDoubleOpMoreEq(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data10,26.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: AthenaStringTypePatternDoubleOpLess
--connection: AthenaApiTests
--input: string some_number = "<20" {pattern: double}
--test: ApiTests:expectTable(AthenaStringTypePatternDoubleOpLess(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data5.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: AthenaStringTypePatternDoubleOpLessEq
--connection: AthenaApiTests
--input: string some_number = "<=20" {pattern: double}
--test: ApiTests:expectTable(AthenaStringTypePatternDoubleOpLessEq(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data5.d42'))
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--CHOICES - should be used for end-to-end tests

--name: AthenaByStringChoices
--input: string country = "France" {choices: ["France", "China", "USA", "Finland"]}
SELECT * FROM mock_data WHERE country = @country;
--end

--name: AthenaByStringChoices
--input: string country = "France" {choices: Query("SELECT DISTINCT country FROM mock_data")}
SELECT * FROM mock_data WHERE country = @country;
--end

--STRING PATTERN

--name: AthenaStringTypePatternStringOpContains
--connection: AthenaApiTests
--input: string first_name = "contains Z" {pattern: string}
--test: ApiTests:expectTable(AthenaStringTypePatternStringOpContains(first_name = 'contains Z'), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data25.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: AthenaStringTypePatternStringOpStartsWith
--connection: AthenaApiTests
--input: string first_name = "starts with W" {pattern: string}
--test: ApiTests:expectTable(AthenaStringTypePatternStringOpStartsWith(first_name='starts with W'), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data23.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: AthenaStringTypePatternStringOpEndsWith
--connection: AthenaApiTests
--input: string first_name = "ends with y" {pattern: string}
--test: ApiTests:expectTable(AthenaStringTypePatternStringOpEndsWith(first_name = 'ends with y'), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data6,23,25.d42'))
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: AthenaStringTypePatternStringOpIn
--connection: AthenaApiTests
--input: string country = "in (Poland, Brazil)" {pattern: string}
--test: ApiTests:expectTable(AthenaStringTypePatternStringOpIn(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data2,5,20.d42'))
SELECT * FROM mock_data WHERE @country(country);
--end

--name: AthenaStringTypePatternStringOpRegex
--connection: AthenaApiTests
--input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
--test: ApiTests:expectTable(AthenaStringTypePatternStringOpRegex(email = 'regex ^([A-Za-z0-9_]+@google.com.au)$'), OpenFile('System:AppData/ApiTests/datasets/tests/athena/data9.d42'))
SELECT * FROM mock_data WHERE @email(email);
--end

--name: AthenaPatternsAllParams
--connection: AthenaApiTests
--input: string first_name = "starts with p" {pattern: string}
--input: string id = ">1" {pattern :int}
--input: bool bool = false
--input: string email = "contains com" {pattern: string}
--input: string some_number = ">20" {pattern: double}
--input: string country = "in (Indonesia)" {pattern: string}
--input: string date = "before 1/1/2022" {pattern: datetime}
--test: ApiTests:expectTable(AthenaPatternsAllParams(first_name = "starts with p", id = ">1", false, email = "contains com", some_number = ">20", country = "in (Indonesia)", date = "before 1/1/2022"), OpenFile("System:AppData/ApiTests/datasets/tests/athena/data13.d42"))
SELECT * FROM mock_data
WHERE @first_name(first_name)
  AND @id(id)
           AND bool = @bool
           AND @email(email)
           AND @some_number(some_number)
           AND @country(country)
           AND @date(date);
--end
