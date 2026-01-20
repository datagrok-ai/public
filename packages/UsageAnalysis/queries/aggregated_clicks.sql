--name: GetAggregatedClicks
--friendlyName: Aggregated Clicks
--connection: System:Datagrok
--meta.cache: true
--meta.cache.invalidateOn: 0 0 * * *

SELECT
    e.description,
    COUNT(*) as count
FROM events e
WHERE e.event_type_id IN (
    SELECT id
    FROM event_types
    WHERE source = 'usage' AND friendly_name = 'click'
)
GROUP BY e.description
ORDER BY count DESC
