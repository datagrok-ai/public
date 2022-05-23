--name: measurementByConceptId
--connection: Mdolotov:PPSCDM_1
--input: int measurement_conc_id
select person_id, value_as_number as measurement_value, measurement_datetime, 
           ROW_NUMBER() OVER(PARTITION BY person_id ORDER BY person_id, measurement_datetime) AS counter
from cdm_pps_prostate_cancer_v2038.measurement
where measurement_concept_id = @measurement_conc_id
limit 100000;