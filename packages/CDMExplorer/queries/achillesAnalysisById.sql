--name: achillesAnalysisByIds
--connection: Demo:test_queries:PostgresCDM
--input: int analysis_id
select * from cdm_synthea_results.achilles_results
where analysis_id = @analysis_id;
