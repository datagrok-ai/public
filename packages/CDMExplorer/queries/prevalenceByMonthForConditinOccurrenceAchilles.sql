--name: PrevalenceByMonthForConditionOccurrenceAchilles
--connection: Demo:test_queries:PostgresCDM
--input: string conceptId
SELECT 
  concept_id                                                AS concept_id,
  concept_name                                              AS concept_name,
 stratum_2                                                AS x_calendar_month,
  -- calendar year, note, there could be blanks
  round(1000 * (1.0 * count_value_num / count_value_denom), 5) AS y_prevalence_1000_pp --prevalence, per 1000 persons
FROM
( select * from 
  (SELECT
     CAST(stratum_1 AS BIGINT) stratum_1_num,
     CAST(stratum_2 AS INT) stratum_2,
     count_value count_value_num
   FROM
     cdm_synthea_results_results WHERE analysis_id = 402 GROUP BY analysis_id, stratum_1, stratum_2, count_value) num
  INNER JOIN
  (SELECT
     CAST(stratum_1 AS INT) stratum_1_denom,
     count_value count_value_denom
   FROM
     cdm_synthea_results_results WHERE analysis_id = 117 GROUP BY analysis_id, stratum_1, count_value) denom
    ON num.stratum_2 = denom.stratum_1_denom
  --calendar year
  INNER JOIN
  cdm_synthea_1.concept c1
ON num.stratum_1_num = c1.concept_id
limit 1000000000) t
WHERE t.concept_id = 200962
ORDER BY CAST(t.stratum_2 AS INT)
