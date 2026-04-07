--name: GetAggregatedClicks
--friendlyName: Aggregated Clicks
--connection: System:Datagrok
--meta.cache: true
--meta.cache.invalidateOn: 0 0 * * *

SELECT
    et.friendly_name as event_type,
    e.description,
    COUNT(*) as count
FROM events e
JOIN event_types et ON e.event_type_id = et.id
WHERE et.source = 'usage' AND et.friendly_name IN ('click', 'menu click', 'dialog show', 'dialog close', 'dialog ok', 'input', 'command', 'navigate')
GROUP BY et.friendly_name, e.description
ORDER BY count DESC
