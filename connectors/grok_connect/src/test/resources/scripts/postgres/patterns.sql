--name: PostgresqlAll
select * from mock_data;
--end



-- INT PATTERN


--name: PostgresqlByInt
--input: int id = 20
SELECT * FROM mock_data WHERE id = @id;
--end

--name: PostgresqlByStringPatternInt
--input: string id = ">28" {pattern: int}
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlByStringPatternInt
--input: string id = ">=29" {pattern: int}
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlByStringPatternInt
--input: string id = "<=1" {pattern: int}
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlByStringPatternInt
--input: string id = "<2" {pattern: int}
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlByStringPatternInt
--input: string id = "in(29, 30)" {pattern: int}
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlByStringPatternInt
--input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlByStringPatternInt
--input: string id = "min-max 29-30" {pattern: int}
SELECT * FROM mock_data WHERE @id(id)
--end

--name: PostgresqlByStringPatternInt
--input: string id = "!=1" {pattern: int}
SELECT * FROM mock_data WHERE @id(id)
--end


--DOUBLE PATTERN

--name: PostgresqlByDouble
--input: double some_number = 510.32
SELECT * FROM mock_data WHERE some_number = @some_number;
--end

--name: PostgresqlByDouble
--input: string some_number = ">975" {pattern: double}
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: PostgresqlByDouble
--input: string some_number = ">=975" {pattern: double}
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: PostgresqlByDouble
--input: string some_number = "<20" {pattern: double}
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: PostgresqlByDouble
--input: string some_number = "<=20" {pattern: double}
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--CHOICES - should be used for end-to-end tests

--name: PostgresqlByStringChoices
--input: string country = 'France' {choices: ["France", "China", "USA", "Finland"]}
SELECT * FROM mock_data WHERE country = @country;
--end

--name: PostgresqlByStringChoices
--input: string country = 'France' {choices: Query("SELECT DISTINCT country FROM mock_data")}
SELECT * FROM mock_data WHERE country = @country;
--end

--STRING PATTERN

--name: PostgresqlByStringPatternString
--input: string first_name = 'contains Z' {pattern: string}
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: PostgresqlByStringPatternString
--input: string first_name = 'starts with W' {pattern: string}
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: PostgresqlByStringPatternString
--input: string first_name = 'ends with y' {pattern: string}
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: PostgresqlByStringPatternString
--input: string country = 'in (Poland, Brazil)' {pattern: string}
SELECT * FROM mock_data WHERE @country(country);
--end

--name: PostgresqlByStringPatternString
--input: string email = 'regex ^([A-Za-z0-9_]+@google.com.au)$' {pattern: string}
SELECT * FROM mock_data WHERE @email(email);
--end

--DATE

--name: PostgresqlByDatetime
--connection: PostgreSQLTest
--input: datetime date = "1/27/2003"
SELECT * FROM mock_data WHERE date >= @date;
--end

--name: PostgresqlByStringPatternDatetime
--input: string date = "today" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: PostgresqlByStringPatternDatetime
--input: string date = "this week" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: PostgresqlByStringPatternDatetime
--input: string date = "this month" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: PostgresqlByStringPatternDatetime
--input: string date = "this year" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: PostgresqlByStringPatternDatetime
--input: string date = "yesterday" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: PostgresqlByStringPatternDatetime
--input: string date = "last year" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: PostgresqlByStringPatternDatetime
--input: string date = "anytime" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: PostgresqlByStringPatternDatetime
--input: string date = "2021-2023" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: PostgresqlByStringPatternDatetime
--input: string date = "before 1/1/2022" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: PostgresqlByStringPatternDatetime
--input: string date = "after 1/1/2022" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: PostgresqlByStringPatternDatetime
--input: string date = "April 2021" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: PostgresqlAll
--input: string first_name = "starts with p" {pattern: string}
--input: string id = ">1" {pattern :int}
--input: bool bool = false
--input: string email = "contains com" {pattern: string}
--input: string some_number = ">20" {pattern: double}
--input: string country = "in (Indonesia)" {pattern: string}
--input: string date = "before 1/1/2022" {pattern: datetime}
SELECT * FROM mock_data
WHERE @first_name(first_name)
  AND @id(id)
           AND bool = @bool
           AND @email(email)
           AND @some_number(some_number)
           AND @country(country)
           AND @date(date);

--end

--name: PostgresqlListParam
--input: list<string> countries = "USA", "France", "China"
SELECT * FROM mock_data WHERE country = ANY(@countries);
