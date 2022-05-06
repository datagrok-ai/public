--name: PostgresqlAll
--connection: PostgreSQLNorthwind
select * from orders;
--end

--name: PostgresqlByInt
--connection: PostgreSQLNorthwind
--input: int orderid
select * from orders where orderid = @orderid;
--end

--name: PostgresqlByStringPatternInt
--connection: PostgreSQLNorthwind
--input: string shipVia {pattern: int}
SELECT * FROM orders WHERE @shipVia(shipVia);
--end

--name: PostgresqlByDouble
--connection: PostgreSQLNorthwind
--input: double freight
SELECT * FROM orders WHERE freight >= @freight;
--end

--name: PostgresqlByStringChoices
--connection: PostgreSQLNorthwind
--input: string shipCountry {choices: ["France", "Germany", "USA", "Finland"]}
SELECT * FROM orders WHERE shipCountry = @shipCountry;
--end

--name: PostgresqlByStringPatternString
--connection: PostgreSQLNorthwind
--input: string shipCity {pattern: string}
SELECT * FROM orders WHERE @shipCity(shipCity);
--end

--name: PostgresqlByBool
--connection: PostgreSQLNorthwind
--input: bool freightLess100
SELECT * FROM orders WHERE ((freight < 100) OR NOT @freightLess100);
--end

--name: PostgresqlByDatetime
--connection: PostgreSQLNorthwind
--input: datetime requiredDate
SELECT * FROM orders WHERE requiredDate >= @requiredDate;
--end

--name: PostgresqlByStringPatternDatetime
--connection: PostgreSQLNorthwind
--input: string orderDate {pattern: datetime}
SELECT * FROM orders WHERE @orderDate(orderDate);
--end
