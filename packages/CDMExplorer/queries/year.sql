--name: year
--connection: Demo:test_queries:PostgresCDM
SELECT 
    stratum_1 AS year, 
    count_value AS countValue
FROM results_pps_prostate_cancer_v2038.achilles_results 
WHERE analysis_id = 3
order by stratum_1 asc
