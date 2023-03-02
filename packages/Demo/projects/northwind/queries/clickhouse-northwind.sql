--name: states
--connection: ClickHouseNorthwind
select distinct stateabbr from usstates order by stateabbr
--end


--name: countries
--connection: ClickHouseNorthwind
select distinct country from customers order by country
--end


--name: products
--connection: ClickHouseNorthwind
select * from products
--end


--name: employees
--connection: ClickHouseNorthwind
select * from employees
--end


--name: customers
--connection: ClickHouseNorthwind
select * from customers
--end


--name: order details by @quantity, @productName, @country
--connection: ClickHouseNorthwind
--input: int quantity = '40'
--input: string productName = 'Manjimup Dried Apples'
--input: string country { choices: northwind:countries }
select
  order_details.orderid,
  order_details.unitprice,
  order_details.quantity,
  products.productname,
  customers.contactname,
  customers.country,
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
  and customers.country = '@country'
  and quantity = @quantity  -- would be nice to be able to construct query on the client side
--end


--name: customers in @country
--connection: ClickHouseNorthwind
--input: string country
select * from customers where country = @country
--end
