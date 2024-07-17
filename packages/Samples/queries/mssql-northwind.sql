--name: MSSQLAll
--connection: MSSQLNorthwind
--meta.testExpectedRows: 830
select * from orders;
--end

--name: MSSQLByInt
--connection: MSSQLNorthwind
--input: int orderid = 10
--meta.testExpectedRows: 1
select * from orders where orderid = @orderid;
--end

--name: MSSQLByStringPatternInt
--connection: MSSQLNorthwind
--input: string shipVia = '>2' {pattern: int}
--meta.testExpectedRows: 255
SELECT * FROM orders WHERE @shipVia(shipVia);
--end

--name: MSSQLByDouble
--connection: MSSQLNorthwind
--input: double freight = 100.1
--meta.testExpectedRows: 187
SELECT * FROM orders WHERE freight >= @freight;
--end

--name: MSSQLByStringChoices
--connection: MSSQLNorthwind
--input: string shipCountry = 'France' {choices: ["France", "Germany", "USA", "Finland"]}
--meta.testExpectedRows: 77
SELECT * FROM orders WHERE shipCountry = @shipCountry;
--end

--name: MSSQLByStringPatternString
--connection: MSSQLNorthwind
--input: string shipCity = 'contains ran' {pattern: string}
--meta.testExpectedRows: 33
SELECT * FROM orders WHERE @shipCity(shipCity);
--end

--name: MSSQLByDatetime
--connection: MSSQLNorthwind
--input: string requiredDate = "after 7/1/1997" {pattern: datetime}
SELECT * FROM orders WHERE @requiredDate(requiredDate);
--end

--name: MSSQLByStringPatternDatetime
--connection: MSSQLNorthwind
--input: string orderDate {pattern: datetime}
SELECT * FROM orders WHERE @orderDate(orderDate);
--end

--name: MSSQLOrders
--friendlyName: Orders
--connection: MSSQLNorthwind
--input: int employeeId = 1
--input: string shipVia = "=3" {pattern: int}
--input: string freight = ">115.0" {pattern: double}
--input: string shipCountry = "USA" {choices: Query("SELECT DISTINCT shipCountry FROM Orders")}
--input: string shipCity = "starts with A" {pattern: string}
--input: bool freightLess1000 = true
--input: string requiredDate = "before 1/1/1997" {pattern: datetime}
--input: string orderDate = "after 1/1/1995" {pattern: datetime}

SELECT * FROM Orders where (employeeId = @employeeId)
                       AND @freight(freight)
    AND @shipVia(shipVia)
   	AND ((freight < 1000) OR @freightLess1000 = 0)
    AND shipCountry = @shipCountry
    AND @shipCity(shipCity)
    AND @orderDate(orderDate)
    AND @requiredDate(requiredDate)
--end

--name: MSSQLProducts
--friendlyName: Products
--connection: MSSQLNorthwind
--input: int ProductID = 7
select * from Products where ProductID = @ProductID
--end
