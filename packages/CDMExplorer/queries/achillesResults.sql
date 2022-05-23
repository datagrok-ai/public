--name: achillesResults
--connection: Mdolotov:PPSCDM_1
select * from results_pps_prostate_cancer_v2038.achilles_results t
join results_pps_prostate_cancer_v2038.achilles_analysis t1 on t.analysis_id = t1.analysis_id
limit 100000;