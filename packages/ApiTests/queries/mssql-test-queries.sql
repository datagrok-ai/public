--name: MSSQLAll
--connection: MSSQLTest
--meta.testExpectedRows: 830
select * from orders;
--end

--name: MSSQLByInt
--connection: MSSQLTest
--input: int orderid = 10330
--meta.testExpectedRows: 1
select * from orders where orderid = @orderid;
--end

--name: MSSQLByStringPatternInt
--connection: MSSQLTest
--input: string shipVia = '>2' {pattern: int}
--meta.testExpectedRows: 255
SELECT * FROM orders WHERE @shipVia(shipVia);
--end

--name: MSSQLByDouble
--connection: MSSQLTest
--input: double freight = 100.1
--meta.testExpectedRows: 187
SELECT * FROM orders WHERE freight >= @freight;
--end

--name: MSSQLByStringChoices
--connection: MSSQLTest
--input: string shipCountry = 'France' {choices: ["France", "Germany", "USA", "Finland"]}
--meta.testExpectedRows: 77
SELECT * FROM orders WHERE shipCountry = @shipCountry;
--end

--name: MSSQLByStringPatternString
--connection: MSSQLTest
--input: string shipCity = 'contains ran' {pattern: string}
--meta.testExpectedRows: 33
SELECT * FROM orders WHERE @shipCity(shipCity);
--end

--name: MSSQLByDatetime
--connection: MSSQLTest
--input: datetime requiredDate
SELECT * FROM orders WHERE requiredDate >= @requiredDate;
--end

--name: MSSQLByStringPatternDatetime
--connection: MSSQLTest
--input: string orderDate {pattern: datetime}
SELECT * FROM orders WHERE @orderDate(orderDate);
--end
