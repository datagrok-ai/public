
--name: acquireCampaignLock
--description: Acquire a lock for a campaign. If the lock is expired, it will be removed and a new lock will be created.
--connection: HitTriage:hitdesign
--input: string appName
--input: string campaignId
--input: string lockedBy {optional: true; nullable: true;}
WITH cleaned AS (
    DELETE FROM hitdesign.campaign_locks 
    WHERE app_name = @appName 
    AND campaign_id = @campaignId 
    AND expires_at < NOW()
)
INSERT INTO  hitdesign.campaign_locks (app_name, campaign_id, expires_at, locked_by)
VALUES (@appName, @campaignId, NOW() + INTERVAL '30 seconds', @lockedBy)
ON CONFLICT (app_name, campaign_id) DO NOTHING
RETURNING app_name, campaign_id, expires_at;

--end

--name: releaseCampaignLock
--description: Release a lock for a campaign.
--connection: HitTriage:hitdesign
--input: string appName
--input: string campaignId
WITH lock_deleted AS (
    DELETE FROM hitdesign.campaign_locks 
    WHERE app_name = @appName
    AND campaign_id = @campaignId
    RETURNING app_name, campaign_id
)
INSERT INTO hitdesign.update_logs (app_name, campaign_id, updated_at)
SELECT app_name, campaign_id, CURRENT_TIMESTAMP
FROM lock_deleted
ON CONFLICT (app_name, campaign_id) 
DO UPDATE SET updated_at = CURRENT_TIMESTAMP
RETURNING updated_at;
--end

--name: getLastModified
--description: Get the last modified time of a campaign lock release (same as save).
--connection: HitTriage:hitdesign
--input: string appName
--input: string campaignId
SELECT updated_at
FROM hitdesign.update_logs
WHERE app_name = @appName
AND campaign_id = @campaignId;
--end

--name: getLastModifiedAll
--description: Get the last save time of every campaign of an app.
--connection: HitTriage:hitdesign
--input: string appName
SELECT campaign_id, updated_at
FROM hitdesign.update_logs
WHERE app_name = @appName;
--end

--name: dropCampaignLock
--description: Drop the caller's lock without logging a save: the campaign was not written under it. A lock another session holds stays.
--connection: HitTriage:hitdesign
--input: string appName
--input: string campaignId
--input: string lockedBy
DELETE FROM hitdesign.campaign_locks
WHERE app_name = @appName
AND campaign_id = @campaignId
AND locked_by = @lockedBy
RETURNING campaign_id;
--end

--name: deleteCampaignRecords
--description: Remove what the database keeps of a deleted campaign: its last save, its lock and its V-iD links.
--connection: HitTriage:hitdesign
--input: string appName
--input: string campaignId
WITH logs AS (
    DELETE FROM hitdesign.update_logs WHERE app_name = @appName AND campaign_id = @campaignId RETURNING 1
), locks AS (
    DELETE FROM hitdesign.campaign_locks WHERE app_name = @appName AND campaign_id = @campaignId RETURNING 1
), vids AS (
    DELETE FROM hitdesign.campaign_vids WHERE app_name = @appName AND campaign_id = @campaignId RETURNING 1
)
SELECT (SELECT COUNT(*) FROM logs) + (SELECT COUNT(*) FROM locks) + (SELECT COUNT(*) FROM vids) AS removed;
--end