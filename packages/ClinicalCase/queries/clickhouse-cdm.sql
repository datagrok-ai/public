--name: First 1000 concepts and domains
--connection: ClickHouseCDM

select
  concept.concept_name,
  concept.domain_id
from
  cdm_synthea.concept
limit 1000

--end
