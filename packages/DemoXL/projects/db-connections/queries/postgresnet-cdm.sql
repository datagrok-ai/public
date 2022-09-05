--name: PostgresNetCDMTest
--friendlyName: CDM
--connection: PostgresNetCDM

select 
  location.city,
  location.state
from
  cdm_synthea_1.location
  
--end