--name: conceptsByIds
--connection: Demo:test_queries:PostgresCDM
--input: string ids = in(0) {pattern: int}
select concept_id, concept_name from cdm_synthea_1.concept
where @ids(concept_id);
