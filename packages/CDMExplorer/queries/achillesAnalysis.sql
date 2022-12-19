--name: achillesAnalysis
--connection: Demo:test_queries:PostgresCDM
select * from cdm_synthea_results.achilles_analysis 
where analysis_id in (select distinct analysis_id from cdm_synthea_results.achilles_results);
