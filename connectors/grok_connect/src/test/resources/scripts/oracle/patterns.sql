--name: OraclePatternsAll
SELECT * FROM mock_data;
--end

-- INT PATTERN

--name: OracleIntTypePatternNone
--input: int id = 20
SELECT * FROM mock_data WHERE id = @id;
--end

--name: OracleStringTypeIntPatternOpMore
--input: string id = ">28" {pattern: int}
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpMoreEq
--input: string id = ">=29" {pattern: int}
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpLessEq
--input: string id = "<=1" {pattern: int}
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpLess
--input: string id = "<2" {pattern: int}
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpIn
--input: string id = "in(29, 30)" {pattern: int}
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpNotIn
--input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpMinMax
--input: string id = "min-max 29-30" {pattern: int}
SELECT * FROM mock_data WHERE @id(id)
--end

--name: OracleStringTypeIntPatternOpNotEq
--input: string id = "!=1" {pattern: int}
SELECT * FROM mock_data WHERE @id(id)
--end


--DOUBLE PATTERN

--name: OracleDoubleTypePatternNone
--input: double some_number = 510.32
SELECT * FROM mock_data WHERE some_number = @some_number;
--end

--name: OracleStringTypePatternDoubleOpMore
--input: string some_number = ">975" {pattern: double}
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: OracleStringTypePatternDoubleOpMoreEq
--input: string some_number = ">=975" {pattern: double}
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: OracleStringTypePatternDoubleOpLess
--input: string some_number = "<20" {pattern: double}
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--name: OracleStringTypePatternDoubleOpLessEq
--input: string some_number = "<=20" {pattern: double}
SELECT * FROM mock_data WHERE @some_number(some_number);
--end

--CHOICES - should be used for end-to-end tests

--name: OracleByStringChoices
--input: string country = "France" {choices: ["France", "China", "USA", "Finland"]}
SELECT * FROM mock_data WHERE country = @country;
--end

--name: OracleByStringChoicesQuery
--input: string country = "France" {choices: Query("SELECT DISTINCT country FROM mock_data")}
SELECT * FROM mock_data WHERE country = @country;
--end

--STRING PATTERN

--name: OracleStringTypePatternStringOpContains
--input: string first_name = "contains Z" {pattern: string}
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: OracleStringTypePatternStringOpStartsWith
--input: string first_name = "starts with W" {pattern: string}
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: OracleStringTypePatternStringOpEndsWith
--input: string first_name = "ends with y" {pattern: string}
SELECT * FROM mock_data WHERE @first_name(first_name);
--end

--name: OracleStringTypePatternStringOpIn
--input: string country = "in (Poland, Brazil)" {pattern: string}
SELECT * FROM mock_data WHERE @country(country);
--end

--name: OracleStringTypePatternStringOpRegex
--input: string email = "regex ^([A-Za-z0-9_]+@google.com.au)$" {pattern: string}
SELECT * FROM mock_data WHERE @email(email);
--end

--DATE

--name: OracleByStringPatternDatetimeOpToday
--input: string date = "today" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: OracleByStringPatternDatetimeOpThisWeek
--input: string date = "this week" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: OracleByStringPatternDatetimeOpThisMonth
--input: string date = "this month" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: OracleByStringPatternDatetimeOpThisYear
--input: string date = "this year" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: OracleByStringPatternDatetimeOpYesterday
--input: string date = "yesterday" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: OracleByStringPatternDatetimeOpLastYear
--input: string date = "last year" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: OracleByStringPatternDatetimeOpAnyTime
--input: string date = "anytime" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: OracleByStringPatternDatetimeRange
--input: string date = "2021-2023" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: OracleByStringPatternDatetimeBefore
--input: string date = "before 1/1/2022" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: OracleByStringPatternDatetimeAfter
--input: string date = "after 1/1/2022" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: OracleByStringPatternDatetime
--input: string date = "April 2021" {pattern: datetime}
SELECT * FROM dates_patterns WHERE @date(date);
--end

--name: OraclePatternsAllParams
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
