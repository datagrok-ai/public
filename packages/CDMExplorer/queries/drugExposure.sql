--name: drugExposure
--connection: Demo:test_queries:PostgresCDM
select t1.person_id, t1.drug_exposure_start_date, t1.drug_exposure_end_date, t2.concept_name as drug
from cdm_synthea_1.drug_exposure t1
join cdm_synthea_1.concept t2
on t1.drug_concept_id = t2.concept_id
limit 10000;
