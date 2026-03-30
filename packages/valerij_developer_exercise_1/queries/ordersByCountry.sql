--name: ordersByCountry
--connection: valerij
--input: string country

select customerid, sum(freight)
from public.orders
where shipcountry = @country
group by customerid