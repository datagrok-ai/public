--name: conditionOccurrence
--connection: Demo:test_queries:PostgresCDM
select t1.person_id, t1.condition_start_date, t1.condition_end_date, t2.concept_name as condition
from cdm_synthea_1.condition_occurrence t1
join cdm_synthea_1.concept t2
on t1.condition_concept_id = t2.concept_id
limit 10000;
