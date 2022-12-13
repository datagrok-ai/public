--name: GenderEthnicityRace
--connection: Demo:test_queries:PostgresCDM
SELECT
  ar1.analysis_id AS analysisId,
  c1.concept_id   AS conceptId,
  c1.concept_name AS conceptName,
  ar1.count_value AS countValue
FROM results_pps_prostate_cancer_v2038.achilles_results ar1
INNER JOIN
cdm_pps_prostate_cancer_v2038.concept c1
ON ar1.stratum_1 = CAST(c1.concept_id AS VARCHAR(255))
WHERE (ar1.analysis_id = 2 AND c1.concept_id IN (8507, 8532)) OR ar1.analysis_id IN (4, 5);
