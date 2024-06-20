--name: ordersByCountry
--connection: agolovko
--input: string country
select customerid, sum(freight)
from public.orders
where shipcountry = @country
group by customerid
-- end
