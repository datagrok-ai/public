--name: getMolecularLiabilityBrowser
--connection: MolecularLiabilityBrowserData:MLB
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
--connection: MolecularLiabilityBrowserData:MLB
select v_id
from public.ptm_observed_v2
--end

--name: getJsonObsByVid
--connection: MolecularLiabilityBrowserData:MLB
--input: string vid
select json
from public.ptm_observed_v2
where v_id = @vid
--end

--name: getJsonComplementByVid
--connection: MolecularLiabilityBrowserData:MLB
--description: PDB positions map to real position  &
--input: string vid = "VR000030945"
select json
from public.jsons_new
where v_id = @vid
--end

--name: getIdMappings
--connection: MolecularLiabilityBrowserData:MLB
select *
from db_v2.gdb_id_mappings
--end

--name: getJsonByVid
--connection: MolecularLiabilityBrowserData:MLB
--description: Chains seq & PTM predictions & Paratop predictions & CDR ranges
--input: string vid = "VR000030945"
select encode(json_data, 'escape')
from db_v2.json_files
where v_id = @vid
--end

--name: getAllJsonsByVid
--connection: MolecularLiabilityBrowserData:MLB
--input: list vids
select *
from db_v2.json_files
where v_id IN @vids
--end

--name: getPdbByVid
--connection: MolecularLiabilityBrowserData:MLB
--input: string vid = "VR000030945"
select encode(pdb_data, 'escape')
from db_v2.pdb_files
where v_id = @vid
--end

--name: getScores
--connection: MolecularLiabilityBrowserData:MLB
--create type core_type_h3l as (h3_length float);
--create type core_type_cdrl as (cdr_length float);
--create type core_type_cdrh as (surface_cdr_hydrophobicity float);
--create type core_type_pos as (positive_cdr_charge float);
--create type core_type_neg as (negative_cdr_charge float);
--create type core_type_sf as (SFvCSP float);
SELECT DISTINCT
ON (v_id) v_id, NGL, gdb_id_mappings, cdr_length, surface_cdr_hydrophobicity, positive_cdr_charge, negative_cdr_charge, SFvCSP
FROM
    (select * FROM db_v2.gdb_id_mappings_v2 t1
    INNER JOIN
    (SELECT v_id as NGL,
    (json_populate_record(null ::core_type_h3l, rightjson::json)).*,
    (json_populate_record(null ::core_type_cdrl, rightjson::json)).*,
    (json_populate_record(null ::core_type_cdrh, rightjson::json)).*,
    (json_populate_record(null ::core_type_pos, rightjson::json)).*,
    (json_populate_record(null ::core_type_neg, rightjson::json)).*,
    (json_populate_record(null ::core_type_sf, rightjson::json)).*
    FROM (WITH t AS (SELECT *, encode(json_tap_data, 'escape') AS myjson FROM db_v2.json_tap_score_files)
    SELECT *, REPLACE(t.myjson, 'SFvCSP', 'sfvcsp') AS rightjson FROM t) t3) t2
    ON t1.v_id=t2.NGL) t4
--end

--name: getVids
--connection: MolecularLiabilityBrowserData:MLB
SELECT jf.v_id
FROM db_v2.pdb_files pf
         JOIN db_v2.json_files jf ON pf.v_id = jf.v_id
--end

--name: listAntigens
--connection: MolecularLiabilityBrowserData:MLB
SELECT ag.id,
       ag.antigen,
       ag.antigen_ncbi_id,
       ag.antigen_gene_symbol
FROM db_v2.antigen as ag
--end

--name: listSchemes
--connection: MolecularLiabilityBrowserData:MLB
SELECT
    s."scheme_id",
    s."scheme"
FROM db_v2.scheme as s
--end

--name: listCdrs
--connection: MolecularLiabilityBrowserData:MLB
SELECT
    cdr."cdr_id",
    cdr."cdr"
FROM db_v2.cdr as cdr
--end

--name: getLayoutBySchemeCdr
--connection: MolecularLiabilityBrowserData:MLB
--input: string scheme
--input: string cdr
SELECT
    sl."region_id",
    sl."scheme",
    sl."cdr",
    sl."type",
    sl."name",
    sl."chain",
    sl."order",
    sl."position_start_name",
    sl."position_end_name"
FROM db_v2.scheme_layout as sl
WHERE sl."scheme" = @scheme and sl."cdr" = @cdr
ORDER BY sl."order", sl."chain"
--end

--name: getMlbByAntigen
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT mlb.v_id,
       mlb.gdb_id_mappings,
       mlb.cdr_length,
       mlb.surface_cdr_hydrophobicity,
       mlb.positive_cdr_charge,
       mlb.negative_cdr_charge,
       mlb."SFvCSP"
FROM public.mlb_main as mlb
         JOIN db_v2.antibody2antigen as ab2ag ON ab2ag.v_id = mlb.v_id
         JOIN db_v2.antigen as ag ON ag.antigen = ab2ag.antigen
WHERE ag.antigen = @antigen;
--end

--name: getTreeByAntigen
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT tr.*
FROM db_v2.tree as tr
         JOIN db_v2.antibody2tree as ab2tr ON ab2tr."CLONE" = tr."CLONE"
         JOIN db_v2.antibody2antigen as ab2ag ON ab2ag."v_id" = ab2tr."v_id"
WHERE ab2ag.antigen = @antigen;
--end

--name: getAnarciChothiaHeavy
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT ma.*
FROM db_v2.antibody_anarci_chothia_heavy as ma
         JOIN db_v2.antibody2antigen as ab2ag ON ab2ag."v_id" = ma."Id"
WHERE ab2ag.antigen = @antigen;
--end

--name: getAnarciChothiaLight
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT ma.*
FROM db_v2.antibody_anarci_chothia_light as ma
         JOIN db_v2.antibody2antigen as ab2ag ON ab2ag."v_id" = ma."Id"
WHERE ab2ag.antigen = @antigen;
--end


--name: getAnarciImgtHeavy
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT ma.*
FROM db_v2.antibody_anarci_imgt_heavy as ma
         JOIN db_v2.antibody2antigen as ab2ag ON ab2ag."v_id" = ma."Id"
WHERE ab2ag.antigen = @antigen;
--end

--name: getAnarciImgtLight
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT ma.*
FROM db_v2.antibody_anarci_imgt_light as ma
         JOIN db_v2.antibody2antigen as ab2ag ON ab2ag."v_id" = ma."Id"
WHERE ab2ag.antigen = @antigen;
--end


--name: getAnarciKabatHeavy
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT ma.*
FROM db_v2.antibody_anarci_kabat_heavy as ma
         JOIN db_v2.antibody2antigen as ab2ag ON ab2ag."v_id" = ma."Id"
WHERE ab2ag.antigen = @antigen;
--end

--name: getAnarciKabatLight
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT ma.*
FROM db_v2.antibody_anarci_kabat_light as ma
         JOIN db_v2.antibody2antigen as ab2ag ON ab2ag."v_id" = ma."Id"
WHERE ab2ag.antigen = @antigen;
--end


--name: getAnarciAhoHeavy
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT ma.*
FROM db_v2.antibody_anarci_aho_heavy as ma
         JOIN db_v2.antibody2antigen as ab2ag ON ab2ag."v_id" = ma."Id"
WHERE ab2ag.antigen = @antigen;
--end

--name: getAnarciAhoLight
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT ma.*
FROM db_v2.antibody_anarci_aho_light as ma
         JOIN db_v2.antibody2antigen as ab2ag ON ab2ag."v_id" = ma."Id"
WHERE ab2ag.antigen = @antigen;
--end
