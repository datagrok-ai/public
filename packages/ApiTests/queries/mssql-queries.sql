--name: MSSQLOrders
--friendlyName: Orders
--connection: MSSQLNorthwind
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
--connection: MSSQLNorthwind
--input: int ProductID = 7
select * from Products where ProductID = @ProductID
--end
