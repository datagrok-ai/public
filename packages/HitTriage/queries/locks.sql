
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