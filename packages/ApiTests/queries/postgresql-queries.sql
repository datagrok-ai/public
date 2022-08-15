--name: PostgresqlOrders
--friendlyName: Orders
--connection: PostgreSQLNorthwind
--tags: unit-test
--input: int employeeId = 5
--input: string shipVia = 3 {pattern: int}
--input: double freight = 10.0
--input: string shipCountry = "France" {choices: Query("SELECT DISTINCT shipCountry FROM Orders")}
--input: string shipCity = "contains r" {pattern: string}
--input: bool freightLess1000 = true
--input: datetime requiredDate = "1/1/1995"
--input: string orderDate = "after 1/1/1995" {pattern: datetime}

SELECT * FROM Orders WHERE (employeeId = @employeeId)
    AND (freight >= @freight)
    AND @shipVia(shipVia)
    AND ((freight < 1000) OR NOT @freightLess1000)
    AND (shipCountry = @shipCountry)
    AND @shipCity(shipCity)
    AND @orderDate(orderDate)
    AND (requiredDate >= @requiredDate)
--end

--name: PostgresqlProducts
--friendlyName: Products
--connection: PostgreSQLNorthwind
--meta.testExpectedRows: 77
select * from Products
--end
