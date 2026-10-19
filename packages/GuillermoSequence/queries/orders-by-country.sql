--name: ordersByCountry
--connection: guillermo
--input: string country = "USA"
select customerid, sum(freight)
from public.orders
where shipcountry = @country
group by customerid
--end
