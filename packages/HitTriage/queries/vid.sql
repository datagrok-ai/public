--name: addMolecule
--friendlyName: Add Molecule
--description: Add a molecule to the dictionary. Returns VID for new or existing molecule.
--connection: HitTriage:hitdesign
--input: string canonicalSmiles [Canonical SMILES of the molecule to register]
--input: string appName [Application the campaign belongs to]
--input: string campaignId [Campaign the molecule is added to]
--input: string createdBy [User registering the molecule]
WITH existing AS (
    SELECT vid FROM hitdesign.vid_dictionary WHERE mh_string = @canonicalSmiles
),
inserted AS (
    INSERT INTO hitdesign.vid_dictionary (mh_string)
    SELECT @canonicalSmiles
    WHERE NOT EXISTS (SELECT 1 FROM existing)
    RETURNING vid
),
all_vids AS (
    SELECT vid FROM existing UNION ALL SELECT vid FROM inserted
),
campaign_insert AS (
    INSERT INTO hitdesign.campaign_vids (app_name, vid, campaign_id, created_by)
    SELECT @appName, vid, @campaignId, @createdBy FROM all_vids
    ON CONFLICT (app_name, vid, campaign_id) DO NOTHING
)
SELECT vid FROM all_vids;
--end

--name: addMolecules
--friendlyName: Add Molecules
--description: Add multiple molecules to the dictionary. Returns VIDs for new or existing molecules.
--connection: HitTriage:hitdesign
--input: list<string> smiles [Canonical SMILES of the molecules to register]
--input: string appName [Application the campaign belongs to]
--input: string campaignId [Campaign the molecules are added to]
--input: string createdBy [User registering the molecules]
WITH input_order AS (
    SELECT s AS mh_string, ROW_NUMBER() OVER () AS ord
    FROM UNNEST(@smiles) AS s
),
existing AS (
    SELECT m.vid, m.mh_string 
    FROM hitdesign.vid_dictionary m
    JOIN input_order i ON m.mh_string = i.mh_string
),
to_insert AS (
    SELECT i.mh_string 
    FROM input_order i
    LEFT JOIN existing e ON i.mh_string = e.mh_string
    WHERE e.mh_string IS NULL
),
inserted AS (
    INSERT INTO hitdesign.vid_dictionary (mh_string)
    SELECT mh_string FROM to_insert
    RETURNING vid, mh_string
),
all_vids AS (
    SELECT vid, mh_string FROM existing 
    UNION ALL 
    SELECT vid, mh_string FROM inserted
),
campaign_inserts AS (
    INSERT INTO hitdesign.campaign_vids (app_name, vid, campaign_id, created_by)
    SELECT @appName, vid, @campaignId, @createdBy FROM all_vids
    ON CONFLICT (app_name, vid, campaign_id) DO NOTHING
)
SELECT v.vid
FROM input_order i
JOIN all_vids v ON v.mh_string = i.mh_string
ORDER BY i.ord;
--end

--name: getMoleculeByVid
--friendlyName: Get Molecule by VID
--description: Get a molecule's canonical SMILES by its VID.
--connection: HitTriage:hitdesign
--input: string vid [Virtual ID (e.g. V000001) to look up]
SELECT mh_string
FROM hitdesign.vid_dictionary
WHERE vid = @vid;
--end

--name: getCampaignsByVid
--friendlyName: Get Campaigns by VID
--description: Get all campaigns and creators for a given VID.
--connection: HitTriage:hitdesign
--input: string vid [Virtual ID (e.g. V000001) to look up]
--input: string appName {optional: true; nullable: true;} [Restrict results to this application]
SELECT campaign_id, created_by, app_name
FROM hitdesign.campaign_vids
WHERE 
    (@appName is NULL or app_name = @appName)
    AND vid = @vid;
--end
