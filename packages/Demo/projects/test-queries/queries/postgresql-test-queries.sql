--name: PostgresqlAll
--connection: PostgreSQLTest
--meta.testExpectedRows: 830
select * from orders;
--end

--name: PostgresqlByInt
--connection: PostgreSQLTest
--input: int orderid = 10330
--meta.testExpectedRows: 1
select * from orders where orderid = @orderid;
--end

--name: PostgresqlByStringPatternInt
--connection: PostgreSQLTest
--input: string shipVia = '>2' {pattern: int}
--meta.testExpectedRows: 255
SELECT * FROM orders WHERE @shipVia(shipVia);
--end

--name: PostgresqlByDouble
--connection: PostgreSQLTest
--input: double freight = 100.1
--meta.testExpectedRows: 187
SELECT * FROM orders WHERE freight >= @freight;
--end

--name: PostgresqlByStringChoices
--connection: PostgreSQLTest
--input: string shipCountry = 'France' {choices: ["France", "Germany", "USA", "Finland"]}
--meta.testExpectedRows: 77
SELECT * FROM orders WHERE shipCountry = @shipCountry;
--end

--name: PostgresqlByStringPatternString
--connection: PostgreSQLTest
--input: string shipCity = 'contains ran' {pattern: string}
--meta.testExpectedRows: 33
SELECT * FROM orders WHERE @shipCity(shipCity);
--end

--name: PostgresqlByBool
--connection: PostgreSQLTest
--input: bool freightLess100 = true
--meta.testExpectedRows: 643
SELECT * FROM orders WHERE ((freight < 100) OR NOT @freightLess100);
--end

--name: PostgresqlByDatetime
--connection: PostgreSQLTest
--input: datetime requiredDate
SELECT * FROM orders WHERE requiredDate >= @requiredDate;
--end

--name: PostgresqlByStringPatternDatetime
--connection: PostgreSQLTest
--input: string orderDate {pattern: datetime}
SELECT * FROM orders WHERE @orderDate(orderDate);
--end
