--name: states
--connection: PostgresNorthwind
select distinct stateabbr from usstates order by stateabbr
--end


--name: countries
--connection: PostgresNorthwind
select distinct country from customers order by country
--end


--name: products
--connection: PostgresNorthwind
select * from products
--end


--name: employees
--connection: PostgresNorthwind
select * from employees
--end


--name: customers
--connection: PostgresNorthwind
select * from customers
--end


--name: order details by @quantity, @productName, @country
--connection: PostgresNorthwind
--input: int quantity = '40'
--input: string productName = 'Manjimup Dried Apples' {suggestions: northwind:productLookup}
--input: string country { choices: northwind:countries }
select
  order_details.orderid,
  order_details.unitprice,
  order_details.quantity,
  products.productname,
  customers.contactname,
  categories.description
from
  order_details
  left join orders on order_details.orderid = orders.orderid
  left join products on order_details.productid = products.productid
  left join customers on orders.customerid = customers.customerid
  left join employees on orders.employeeid = employees.employeeid
  left join shippers on orders.shipvia = shippers.shipperid
  left join suppliers on products.supplierid = suppliers.supplierid
  left join categories on products.categoryid = categories.categoryid
where
  products.productname like '%@productName%'
  and country = '@country'
  and quantity = @quantity  -- would be nice to be able to construct query on the client side
--end


--name: customers in @country
--connection: PostgresNorthwind
--input: string country
select * from customers where country = @country
--end


--name: productLookup
--connection: PostgresNorthwind
--input: string lookup
select productname from products where productname like (@lookup || '%')


--name: productDetails
--connection: PostgresNorthwind
--input: string productName {suggestions: northwind:productLookup}
select * from products where productname = @productName


--name: OrdersByEmployee
--friendlyName: OrdersByEmployee
--connection: PostgresNorthwind
--input: string shipCountry = "Argentina" {choices: Query("SELECT DISTINCT shipCountry FROM orders")}
--input: string shipCity {choices: Query("SELECT DISTINCT shipcity FROM orders WHERE shipCountry = @shipCountry")}
--input: string customerId {choices: Query("SELECT DISTINCT customerid FROM orders WHERE shipCity = @shipCity")}
--input: string employee {choices: Query("SELECT DISTINCT lastName FROM employees INNER JOIN orders ON employees.employeeId = orders.employeeId WHERE customerId = @customerId")}

SELECT DISTINCT orderid
FROM orders
INNER JOIN employees
ON orders.employeeId = employees.employeeId
WHERE lastName = @employee

--end