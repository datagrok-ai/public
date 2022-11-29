--name: getMolecularLiabilityBrowser
--connection: MolecularLiabilityBrowserData:MLB
--description: Used by old MolecularLiabilityBrowser
SELECT v_id,
       gdb_id_mappings,
       cdr_length,
       surface_cdr_hydrophobicity,
       positive_cdr_charge,
       negative_cdr_charge,
       SFvCSP
FROM mlb.mlb_main
--end

--name: getObservedPtmVids
--connection: MolecularLiabilityBrowserData:MLB
select v_id
from mlb.ptm_observed_v2
--end

--name: getJsonObsByVid
--connection: MolecularLiabilityBrowserData:MLB
--input: string vid
select json
from mlb.ptm_observed_v2
where v_id = @vid
--end

--name: getJsonComplementByVid
--connection: MolecularLiabilityBrowserData:MLB
--description: PDB positions map to real position  &
--input: string vid = "VR000030945"
select json
from mlb.jsons_new
where v_id = @vid
--end

--name: getJsonByVid
--connection: MolecularLiabilityBrowserData:MLB
--description: Chains seq & PTM predictions & Paratop predictions & CDR ranges
--input: string vid = "VR000030945"
select encode(json_data, 'escape')
from db_v2.json_files
where v_id = @vid
--end

--name: getPdbByVid
--connection: MolecularLiabilityBrowserData:MLB
--input: string vid = "VR000030945"
select encode(pdb_data, 'escape')
from db_v2.pdb_files
where v_id = @vid
--end

--name: getListVersion
--connection: MolecularLiabilityBrowserData:MLB
SELECT lv.list_id, lv.name, lv.version
FROM mlb.list_version as lv;
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
       ag.antigen_gene_symbol,
       (SELECT count(*) FROM mlb.tree2 as tr2 WHERE tr2.antigen = ag.antigen AND tr2."CLONE" <> 'REPERTOIRE') as clones
FROM mlb.antigen as ag
         JOIN mlb.antibody2antigen as ab2ag ON ab2ag.antigen = ag.antigen
GROUP BY ag.id, ag.antigen;
--end

--name: listSchemes
--connection: MolecularLiabilityBrowserData:MLB
SELECT s."scheme_id",
       s."scheme"
FROM mlb.scheme as s
--end

--name: listCdrs
--connection: MolecularLiabilityBrowserData:MLB
SELECT cdr."cdr_id",
       cdr."cdr"
FROM mlb.cdr as cdr
--end

--name: getLayoutBySchemeCdr
--connection: MolecularLiabilityBrowserData:MLB
--input: string scheme
--input: string cdr
SELECT sl."region_id",
       sl."scheme",
       sl."cdr",
       sl."type",
       sl."name",
       sl."chain",
       sl."order",
       sl."position_start_name",
       sl."position_end_name"
FROM mlb.scheme_layout as sl
WHERE sl."scheme" = @scheme
  and sl."cdr" = @cdr
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
FROM mlb.mlb_main as mlb
         JOIN mlb.antibody2antigen as ab2ag ON ab2ag.v_id = mlb.v_id
         JOIN mlb.antigen as ag ON ag.antigen = ab2ag.antigen
WHERE ag.antigen = @antigen;
--end

--name: getTreeByAntigen
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT tr2.*
FROM mlb.tree2 as tr2
WHERE (tr2.antigen = @antigen) AND (tr2."CLONE" <> 'REPERTOIRE');
--end

--name: getAnarciChothiaHeavy
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT ma.*
FROM mlb.antibody_anarci_chothia_heavy as ma
         JOIN mlb.antibody2antigen as ab2ag ON ab2ag."v_id" = ma."Id"
WHERE ab2ag.antigen = @antigen;
--end

--name: getAnarciChothiaLight
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT ma.*
FROM mlb.antibody_anarci_chothia_light as ma
         JOIN mlb.antibody2antigen as ab2ag ON ab2ag."v_id" = ma."Id"
WHERE ab2ag.antigen = @antigen;
--end


--name: getAnarciImgtHeavy
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT ma.*
FROM mlb.antibody_anarci_imgt_heavy as ma
         JOIN mlb.antibody2antigen as ab2ag ON ab2ag."v_id" = ma."Id"
WHERE ab2ag.antigen = @antigen;
--end

--name: getAnarciImgtLight
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT ma.*
FROM mlb.antibody_anarci_imgt_light as ma
         JOIN mlb.antibody2antigen as ab2ag ON ab2ag."v_id" = ma."Id"
WHERE ab2ag.antigen = @antigen;
--end


--name: getAnarciKabatHeavy
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT ma.*
FROM mlb.antibody_anarci_kabat_heavy as ma
         JOIN mlb.antibody2antigen as ab2ag ON ab2ag."v_id" = ma."Id"
WHERE ab2ag.antigen = @antigen;
--end

--name: getAnarciKabatLight
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT ma.*
FROM mlb.antibody_anarci_kabat_light as ma
         JOIN mlb.antibody2antigen as ab2ag ON ab2ag."v_id" = ma."Id"
WHERE ab2ag.antigen = @antigen;
--end


--name: getAnarciAhoHeavy
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT ma.*
FROM mlb.antibody_anarci_aho_heavy as ma
         JOIN mlb.antibody2antigen as ab2ag ON ab2ag."v_id" = ma."Id"
WHERE ab2ag.antigen = @antigen;
--end

--name: getAnarciAhoLight
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
SELECT ma.*
FROM mlb.antibody_anarci_aho_light as ma
         JOIN mlb.antibody2antigen as ab2ag ON ab2ag."v_id" = ma."Id"
WHERE ab2ag.antigen = @antigen;
--end


--name: getPredictedPtmByAntigen
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
select
    ptm_pred.*
from
    mlb.ptm_predicted_v2 as ptm_pred
         JOIN mlb.antibody2antigen as ab2ag ON ab2ag.v_id = ptm_pred.v_id
         JOIN mlb.antigen as ag ON ag.antigen = ab2ag.antigen
WHERE ag.antigen = @antigen;
--end

--name: getObservedPtmByAntigen
--connection: MolecularLiabilityBrowserData:MLB
--input: string antigen
select
    ptm_obs.*
from
    mlb.ptm_observed_v2 as ptm_obs
         JOIN mlb.antibody2antigen as ab2ag ON ab2ag.v_id = ptm_obs.v_id
         JOIN mlb.antigen as ag ON ag.antigen = ab2ag.antigen
WHERE ag.antigen = @antigen;
--end


--name: get3D
--connection: MolecularLiabilityBrowserData:MLB
--input: string vid = "VR000030945"
SELECT
    (select encode(json_data, 'escape')
        from db_v2.json_files
        where v_id = @vid) as "json",
       (select encode(pdb_data, 'escape')
        from db_v2.pdb_files
        where v_id = @vid) as "pdb",
       (select json
        from mlb.jsons_new
        where v_id = @vid) as "real_nums",
       (select json
        from db_v2.json_files_observed
        where v_id = @vid
                              limit 1) as "obs_ptm";
--end
