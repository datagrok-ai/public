--name: conceptsByNumbers
--connection: Mdolotov:PPSCDM_1
--string ids
select concept_id, concept_name from cdm_pps_prostate_cancer_v2038.achilles_results
where concept_id in (@ids);