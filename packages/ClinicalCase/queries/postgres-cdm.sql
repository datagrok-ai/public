--name: Сities and States
--connection: PostgresCDM

select
  location.city,
  location.state
from
  cdm_synthea_1.location

--end
