--name: cohorts
--connection: Demo:test_queries:PostgresCDM
select t1.person_id, t2.cohort_definition_id, t2.cohort_start_date, t2.cohort_end_date, t3.cohort_definition_name
from cdm_synthea_1.person t1
left join cdm_synthea_1.cohort t2 on t1.person_id = t2.subject_id
left join cdm_synthea_1.cohort_definition t3 on t2.cohort_definition_id = t3.cohort_definition_id
