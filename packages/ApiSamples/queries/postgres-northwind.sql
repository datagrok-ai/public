--name: Countries
--connection: PostgresNorthwind
select distinct country from customers order by country
--end

--name: CustomersInCountry
--connection: PostgresNorthwind
--input: string country = "France"
select * from customers where country = @country
--end
