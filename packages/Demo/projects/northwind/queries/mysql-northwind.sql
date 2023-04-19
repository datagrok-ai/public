--name: Orders
--friendlyName: Orders
--connection: MySQLNorthwind
--input: int employeeId = 5
--input: string shipVia = "= 3" {pattern: int}
--input: double freight = 10.0
--input: string shipCountry = "france" {choices: query("select distinct shipCountry from orders")}
--input: string shipCity = "starts with r" {pattern: string}
--input: bool freightLess1000 = true
--input: datetime requiredDate = "1/1/1995"
--input: string orderDate = "after 1/1/1995" {pattern: datetime}

SELECT *
FROM orders
WHERE (employeeId = @employeeId)
  and (freight >= @freight)
  and @shipVia(shipVia)
    and ((freight < 1000) or not @freightLess1000)
    and (shipCountry = @shipCountry)
    and @shipCity(shipCity)
    and @orderDate(orderDate)
    and (requiredDate >= @requiredDate)

--end

--name: Products
--friendlyName: Products
--connection: MySQLNorthwind

SELECT *
FROM products

--end

