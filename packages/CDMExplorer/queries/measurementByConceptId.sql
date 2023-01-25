--name: measurementByConceptId
--connection: Demo:test_queries:PostgresCDM
--input: int measurement_conc_id
select person_id, value_as_number as measurement_value, measurement_datetime, measurement_concept_id,
EXTRACT(YEAR FROM measurement_datetime) m_year, 
EXTRACT(MONTH FROM measurement_datetime) m_month, 
EXTRACT(WEEK FROM measurement_datetime) m_week,
           ROW_NUMBER() OVER(PARTITION BY person_id ORDER BY measurement_concept_id, person_id, measurement_datetime) AS counter
from cdm_synthea_1.measurement
where measurement_concept_id = @measurement_conc_id
limit 100000;
