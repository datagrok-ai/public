--name: MSSQLAll
--connection: MSSQLNorthwind
select * from orders;
--end

--name: MSSQLByInt
--connection: MSSQLNorthwind
--input: int orderid
select * from orders where orderid = @orderid;
--end

--name: MSSQLByStringPatternInt
--connection: MSSQLNorthwind
--input: string shipVia {pattern: int}
SELECT * FROM orders WHERE @shipVia(shipVia);
--end

--name: MSSQLByDouble
--connection: MSSQLNorthwind
--input: double freight
SELECT * FROM orders WHERE freight >= @freight;
--end

--name: MSSQLByStringChoices
--connection: MSSQLNorthwind
--input: string shipCountry {choices: ["France", "Germany", "USA", "Finland"]}
SELECT * FROM orders WHERE shipCountry = @shipCountry;
--end

--name: MSSQLByStringPatternString
--connection: MSSQLNorthwind
--input: string shipCity {pattern: string}
SELECT * FROM orders WHERE @shipCity(shipCity);
--end

--name: MSSQLByDatetime
--connection: MSSQLNorthwind
--input: datetime requiredDate
SELECT * FROM orders WHERE requiredDate >= @requiredDate;
--end

--name: MSSQLByStringPatternDatetime
--connection: MSSQLNorthwind
--input: string orderDate {pattern: datetime}
SELECT * FROM orders WHERE @orderDate(orderDate);
--end