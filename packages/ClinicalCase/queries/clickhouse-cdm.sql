--name: Concepts and Domains
--connection: ClickHouseCDM

select
  concept.concept_name,
  concept.domain_id
from
  cdm_synthea.concept

--end
