
--name: davitQuerry
--connection: NorthwindTest
--input: string country
--output: dataframe result

select customerid, sum(freight)
from public.orders
where shipcountry = @country
group by customerid
--end
