--name: measurementTypes
--connection: Mdolotov:PPSCDM_1
select distinct t1.measurement_concept_id as concept_id, t2.concept_name as measurement
from cdm_pps_prostate_cancer_v2038.measurement t1
join cdm_pps_prostate_cancer_v2038.concept t2
on t1.measurement_concept_id = t2.concept_id