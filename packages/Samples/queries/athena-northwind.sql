--name: Orders
--friendlyName: Orders
--connection: AthenaNorthwind

--input: int employeeId = 5
--input: string shipVia = "= 3" {pattern: int}
--input: double freight = 10.0
--input: string shipCountry = "France" {choices: Query("SELECT DISTINCT shipCountry FROM northwind.orders")}
--input: string shipCity = "starts with r" {pattern: string}
--input: bool freightLess1000 = true
--input: string orderDate = "after 1/1/1995" {pattern: datetime}

SELECT * FROM northwind.orders WHERE (employeeId = @employeeId)
    AND (freight >= @freight)
    AND @shipVia(shipVia)
    AND ((freight < 1000) OR NOT @freightLess1000)
    AND (shipCountry = @shipCountry)
    AND @shipCity(shipCity)
    AND @orderDate(orderDate)

--end