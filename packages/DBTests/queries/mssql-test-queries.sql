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

--name: MSSQLOrders
--friendlyName: Orders
--connection: MSSQLTest
--input: int employeeId = 5
--input: string shipVia = "=3" {pattern: int}
--input: double freight = 10.0
--input: string shipCountry = "France" {choices: Query("SELECT DISTINCT shipCountry FROM Orders")}
--input: string shipCity = "starts with r" {pattern: string}
--input: bool freightLess1000 = true
--input: datetime requiredDate = "1/1/1995"
--input: string orderDate = "after 1/1/1995" {pattern: datetime}

SELECT * FROM Orders WHERE (employeeId = @employeeId)
    AND (freight >= @freight)
    AND @shipVia(shipVia)
    AND ((freight < 1000) OR @freightLess1000 = 0)
    AND (shipCountry = @shipCountry)
    AND @shipCity(shipCity)
    AND @orderDate(orderDate)
    AND (requiredDate >= @requiredDate)
--end

--name: MSSQLProducts
--friendlyName: Products
--connection: MSSQLTest
--input: int ProductID = 7
select * from Products where ProductID = @ProductID
--end
