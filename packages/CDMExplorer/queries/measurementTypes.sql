--name: measurementTypes
--connection: Demo:test_queries:PostgresCDM
select distinct t1.measurement_concept_id as concept_id, t2.concept_name as measurement
from cdm_synthea_1.measurement t1
join cdm_synthea_1.concept t2
on t1.measurement_concept_id = t2.concept_id
