--name: getMolecularLiabilityBrowser
--connection: MLB
SELECT v_id,
       gdb_id_mappings,
       cdr_length,
       surface_cdr_hydrophobicity,
       positive_cdr_charge,
       negative_cdr_charge,
       SFvCSP
FROM public.mlb_main
-- NGL column removed from query 2022-04-08 (atanas/lstolbov)
--end

--name: getObservedPtmVids
--connection: MLB
select v_id
from public.ptm_observed_v2
--end

--name: getJsonObsByVid
--connection: MLB
--input: string vid
select json
from public.ptm_observed_v2
where v_id = @vid
--end

--name: getJsonComplementByVid
--connection: MLB
--input: string vid = "VR000030945"
select json
from public.jsons_new
where v_id = @vid
--end

--name: getIdMappings
--connection: MLB
select *
from db_v2.gdb_id_mappings
--end

--name: getJsonByVid
--connection: MLB
--input: string vid = "VR000030945"
select encode(json_data, 'escape')
from db_v2.json_files
where v_id = @vid
--end

--name: getAllJsonsByVid
--connection: MLB
--input: list vids
select *
from db_v2.json_files
where v_id IN @vids
--end

--name: getPdbByVid
--connection: MLB
--input: string vid = "VR000030945"
select encode(pdb_data, 'escape')
from db_v2.pdb_files
where v_id = @vid
--end

--name: getScores
--connection: MLB
--create type core_type_h3l as (h3_length float);
--create type core_type_cdrl as (cdr_length float);
--create type core_type_cdrh as (surface_cdr_hydrophobicity float);
--create type core_type_pos as (positive_cdr_charge float);
--create type core_type_neg as (negative_cdr_charge float);
--create type core_type_sf as (SFvCSP float);
SELECT DISTINCT ON (v_id) v_id, NGL, gdb_id_mappings, cdr_length, surface_cdr_hydrophobicity, positive_cdr_charge, negative_cdr_charge, SFvCSP  FROM
	(select * FROM db_v2.gdb_id_mappings_v2 t1
  	INNER JOIN
	(SELECT v_id as NGL,
		(json_populate_record(null::core_type_h3l, rightjson::json)).*,
        (json_populate_record(null::core_type_cdrl, rightjson::json)).*,
        (json_populate_record(null::core_type_cdrh, rightjson::json)).*,
        (json_populate_record(null::core_type_pos, rightjson::json)).*,
        (json_populate_record(null::core_type_neg, rightjson::json)).*,
        (json_populate_record(null::core_type_sf, rightjson::json)).*
		FROM (WITH t AS (SELECT *, encode(json_tap_data, 'escape') AS myjson FROM db_v2.json_tap_score_files)
			SELECT *, REPLACE(t.myjson, 'SFvCSP', 'sfvcsp') AS rightjson FROM t) t3) t2
	ON t1.v_id=t2.NGL) t4
--end

--name: getVids
--connection: MLB
SELECT
	jf.v_id
FROM
	db_v2.pdb_files pf
	JOIN db_v2.json_files jf ON pf.v_id = jf.v_id
--end

--name: listAntigens
--connection: MLB
SELECT
    ag.id, ag.antigen, ag.antigen_ncbi_id, ag.antigen_gene_symbol
FROM
    db_v2.antigen as ag
--end

