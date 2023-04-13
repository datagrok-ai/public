--name: Orders
--friendlyName: Orders
--connection: MySQLNorthwind
--input: int employeeid = 5
--input: string shipvia = "= 3" {pattern: int}
--input: double freight = 10.0
--input: string shipcountry = "france" {choices: query("select distinct shipcountry from orders")}
--input: string shipcity = "starts with r" {pattern: string}
--input: bool freightless1000 = true
--input: datetime requireddate = "1/1/1995"
--input: string orderdate = "after 1/1/1995" {pattern: datetime}

SELECT *
FROM orders
WHERE (employeeid = @employeeid)
  and (freight >= @freight)
  and @shipvia(shipvia)
    and ((freight < 1000) or not @freightless1000)
    and (shipcountry = @shipcountry)
    and @shipcity(shipcity)
    and @orderdate(orderdate)
    and (requireddate >= @requireddate)

--end

--name: Products
--friendlyName: Products
--connection: MySQLNorthwind

SELECT *
FROM products

--end

--name: OrdersByEmployeeID
--friendlyName: OrdersByEmployeeID
--connection: MySQLNorthwind
--input: string shipcountry = "Argentina" {choices: Query("SELECT DISTINCT shipcountry FROM orders")}
--input: string shipcity {choices: Query("SELECT DISTINCT shipcity FROM orders WHERE shipcountry = @shipcountry")}
--input: string customerid {choices: Query("SELECT DISTINCT customerid FROM orders WHERE shipcity = @shipcity")}
--input: string employeeid {choices: Query("SELECT DiSTINCT employeeid FROM orders WHERE customerid = @customerid")}

SELECT DISTINCT orderid
FROM orders
WHERE (employeeid = CAST(@employeeid AS INTEGER))

--end