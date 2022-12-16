--name: person
--connection: Demo:test_queries:PostgresCDM
select t1.*, t2.concept_name as gender, t3.concept_name as race, t4.concept_name as ethnicity
from cdm_synthea_1.person t1
left join cdm_synthea_1.concept t2 on t1.gender_concept_id = t2.concept_id
left join cdm_synthea_1.concept t3 on t1.race_concept_id = t3.concept_id
left join cdm_synthea_1.concept t4 on t1.ethnicity_concept_id = t4.concept_id
limit 100000;
