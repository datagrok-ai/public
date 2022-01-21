--name: PostgresNetAll
--connection: PostgresNetNorthwind
select * from orders;
--end

--name: PostgresNetByInt
--connection: PostgresNetNorthwind
--input: int orderid
select * from orders where orderid = @orderid;
--end

--name: PostgresNetByStringPatternInt
--connection: PostgresNetNorthwind
--input: string shipVia {pattern: int}
SELECT * FROM orders WHERE @shipVia(shipVia);
--end

--name: PostgresNetByDouble
--connection: PostgresNetNorthwind
--input: double freight
SELECT * FROM orders WHERE freight >= @freight;
--end

--name: PostgresNetByStringChoices
--connection: PostgresNetNorthwind
--input: string shipCountry {choices: ["France", "Germany", "USA", "Finland"]}
SELECT * FROM orders WHERE shipCountry = @shipCountry;
--end

--name: PostgresNetByStringPatternString
--connection: PostgresNetNorthwind
--input: string shipCity {pattern: string}
SELECT * FROM orders WHERE @shipCity(shipCity);
--end

--name: PostgresNetByBool
--connection: PostgresNetNorthwind
--input: bool freightLess100
SELECT * FROM orders WHERE ((freight < 100) OR NOT @freightLess100);
--end

--name: PostgresNetByDatetime
--connection: PostgresNetNorthwind
--input: datetime requiredDate
SELECT * FROM orders WHERE requiredDate >= @requiredDate;
--end

--name: PostgresNetByStringPatternDatetime
--connection: PostgresNetNorthwind
--input: string orderDate {pattern: datetime}
SELECT * FROM orders WHERE @orderDate(orderDate);
--end