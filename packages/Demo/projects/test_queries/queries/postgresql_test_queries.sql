--name: PostgresqlAll
--connection: PostgresqlNorthwind
select * from orders;
--end

--name: PostgresqlByInt
--connection: PostgresqlNorthwind
--input: int orderid
select * from orders where orderid = @orderid;
--end

--name: PostgresqlByStringPatternInt
--connection: PostgresqlNorthwind
--input: string shipVia {pattern: int}
SELECT * FROM orders WHERE @shipVia(shipVia);
--end

--name: PostgresqlByDouble
--connection: PostgresqlNorthwind
--input: double freight
SELECT * FROM orders WHERE freight >= @freight;
--end

--name: PostgresqlByStringChoices
--connection: PostgresqlNorthwind
--input: string shipCountry {choices: ["France", "Germany", "USA", "Finland"]}
SELECT * FROM orders WHERE shipCountry = @shipCountry;
--end

--name: PostgresqlByStringPatternString
--connection: PostgresqlNorthwind
--input: string shipCity {pattern: string}
SELECT * FROM orders WHERE @shipCity(shipCity);
--end

--name: PostgresqlByBool
--connection: PostgresqlNorthwind
--input: bool freightLess100
SELECT * FROM orders WHERE ((freight < 100) OR NOT @freightLess100);
--end

--name: PostgresqlByDatetime
--connection: PostgresqlNorthwind
--input: datetime requiredDate
SELECT * FROM orders WHERE requiredDate >= @requiredDate;
--end

--name: PostgresqlByStringPatternDatetime
--connection: PostgresqlNorthwind
--input: string orderDate {pattern: datetime}
SELECT * FROM orders WHERE @orderDate(orderDate);
--end