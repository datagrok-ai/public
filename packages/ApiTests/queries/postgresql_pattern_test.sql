--name: PostgresqlPatternsAll
--connection: PostgreSQLApiTests
--test: ApiTests:compareTables(PostgresqlPatternsAll(), ApiTests:getTable("data1-30.csv", path="ApiTests/datasets/tests"))
SELECT * FROM mock_data;
--end

-- INT PATTERN

--name: PostgresqlIntTypePatternNone
--connection: PostgreSQLApiTests
--input: int id = 20
--test: ApiTests:compareTables(PostgresqlIntTypePatternNone(), ApiTests:getTable('data20.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE id = @id;
--end

--name: PostgresqlStringTypeIntPatternOpMore
--connection: PostgreSQLApiTests
--input: string id = ">28" {pattern: int}
--test: ApiTests:compareTables(PostgresqlStringTypeIntPatternOpMore(), ApiTests:getTable('data29-30.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlStringTypeIntPatternOpMoreEq
--connection: PostgreSQLApiTests
--input: string id = ">=29" {pattern: int}
--test: ApiTests:compareTables(PostgresqlStringTypeIntPatternOpMoreEq(), ApiTests:getTable('data29-30.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlStringTypeIntPatternOpLessEq
--connection: PostgreSQLApiTests
--input: string id = "<=1" {pattern: int}
--test: ApiTests:compareTables(PostgresqlStringTypeIntPatternOpLessEq(), ApiTests:getTable('data1.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlStringTypeIntPatternOpLess
--connection: PostgreSQLApiTests
--input: string id = "<2" {pattern: int}
--test: ApiTests:compareTables(PostgresqlStringTypeIntPatternOpLess(), ApiTests:getTable('data1.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlStringTypeIntPatternOpIn
--connection: PostgreSQLApiTests
--input: string id = "in(29, 30)" {pattern: int}
--test: ApiTests:compareTables(PostgresqlStringTypeIntPatternOpIn(), ApiTests:getTable('data29-30.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlStringTypeIntPatternOpNotIn
--connection: PostgreSQLApiTests
--input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
--test: ApiTests:compareTables(PostgresqlStringTypeIntPatternOpNotIn(), ApiTests:getTable('data1-20.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlStringTypeIntPatternOpMinMax
--connection: PostgreSQLApiTests
--input: string id = "min-max 29-30" {pattern: int}
--test: ApiTests:compareTables(PostgresqlStringTypeIntPatternOpMinMax(), ApiTests:getTable('data29-30.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlStringTypeIntPatternOpNotEq
--connection: PostgreSQLApiTests
--input: string id = "!=1" {pattern: int}
--test: ApiTests:compareTables(PostgresqlStringTypeIntPatternOpNotEq(), ApiTests:getTable('data2-30.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE @id(id)
--end


--DOUBLE PATTERN

--name: PostgresqlDoubleTypePatternNone
--connection: PostgreSQLApiTests
--input: double some_number = 510.32
--test: ApiTests:compareTables(PostgresqlDoubleTypePatternNone(), ApiTests:getTable('data1.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE some_number = @some_number;
--end

--name: PostgresqlStringTypePatternDoubleOpMore
--connection: PostgreSQLApiTests
--input: string some_number = ">975" {pattern: double}
--test: ApiTests:compareTables(PostgresqlStringTypePatternDoubleOpMore(), ApiTests:getTable('data10,26.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: PostgresqlStringTypePatternDoubleOpMoreEq
--connection: PostgreSQLApiTests
--input: string some_number = ">=975" {pattern: double}
--test: ApiTests:compareTables(PostgresqlStringTypePatternDoubleOpMoreEq(), ApiTests:getTable('data10,26.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: PostgresqlStringTypePatternDoubleOpLess
--connection: PostgreSQLApiTests
--input: string some_number = "<20" {pattern: double}
--test: ApiTests:compareTables(PostgresqlStringTypePatternDoubleOpLess(), ApiTests:getTable('data5.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: PostgresqlStringTypePatternDoubleOpLessEq
--connection: PostgreSQLApiTests
--input: string some_number = "<=20" {pattern: double}
--test: ApiTests:compareTables(PostgresqlStringTypePatternDoubleOpLessEq(), ApiTests:getTable('data5.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--CHOICES - should be used for end-to-end tests

--name: PostgresqlByStringChoices
--input: string country = "France" {choices: ["France", "China", "USA", "Finland"]}
SELECT * FROM mock_data WHERE country = @country;
--end

--name: PostgresqlByStringChoices
--input: string country = "France" {choices: Query("SELECT DISTINCT country FROM mock_data")}
SELECT * FROM mock_data WHERE country = @country;
--end

--STRING PATTERN

--name: PostgresqlStringTypePatternStringOpContains
--connection: PostgreSQLApiTests
--input: string first_name = "contains Z" {pattern: string}
--test: ApiTests:compareTables(PostgresqlStringTypePatternStringOpContains(), ApiTests:getTable('data25.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: PostgresqlStringTypePatternStringOpStartsWith
--connection: PostgreSQLApiTests
--input: string first_name = "starts with W" {pattern: string}
--test: ApiTests:compareTables(PostgresqlStringTypePatternStringOpStartsWith(), ApiTests:getTable('data23.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: PostgresqlStringTypePatternStringOpEndsWith
--connection: PostgreSQLApiTests
--input: string first_name = "ends with y" {pattern: string}
--test: ApiTests:compareTables(PostgresqlStringTypePatternStringOpEndsWith(), ApiTests:getTable('data6,23,25.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: PostgresqlStringTypePatternStringOpIn
--connection: PostgreSQLApiTests
--input: string country = "in (Poland, Brazil)" {pattern: string}
--test: ApiTests:compareTables(PostgresqlStringTypePatternStringOpIn(), ApiTests:getTable('data2,5,20.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE @country(country);
--end

--name: PostgresqlStringTypePatternStringOpRegex
--connection: PostgreSQLApiTests
--input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
--test: ApiTests:compareTables(PostgresqlStringTypePatternStringOpRegex(), ApiTests:getTable('data9.csv', path='ApiTests/datasets/tests'))
SELECT * FROM mock_data WHERE @email(email);
--end

--name: PostgresqlPatternsAllParams
--connection: PostgreSQLApiTests
--input: string first_name = "starts with p" {pattern: string}
--input: string id = ">1" {pattern :int}
--input: bool bool = false
--input: string email = "contains com" {pattern: string}
--input: string some_number = ">20" {pattern: double}
--input: string country = "in (Indonesia)" {pattern: string}
--input: string date = "before 1/1/2022" {pattern: datetime}
--test: ApiTests:compareTables(PostgresqlPatternsAllParams(), ApiTests:getTable("data13.csv", path="ApiTests/datasets/tests"))
SELECT * FROM mock_data
WHERE @first_name(first_name)
  AND @id(id)
           AND bool = @bool
           AND @email(email)
           AND @some_number(some_number)
           AND @country(country)
           AND @date(date);

--end
