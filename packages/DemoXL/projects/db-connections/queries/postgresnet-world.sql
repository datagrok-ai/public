--name: PostgresNetWorld
--friendlyName: World
--connection: PostgresNetWorld

select 
  column_domain_usage.column_name,
  column_domain_usage.domain_catalog,
  column_domain_usage.domain_name,
  column_domain_usage.table_schema
from
  information_schema.column_domain_usage

--end