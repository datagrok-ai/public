--name: GetAggregatedClicks
--friendlyName: Aggregated Clicks
--connection: System:Datagrok
--input: string date { pattern: datetime }

SELECT
    et.friendly_name as event_type,
    e.description,
    COUNT(*)::int as count
FROM events e
JOIN event_types et ON e.event_type_id = et.id
WHERE et.source = 'usage'
  AND et.friendly_name IN ('click', 'menu click', 'dialog show', 'dialog close', 'dialog ok', 'input', 'command', 'navigate')
  AND @date(e.event_time)
GROUP BY et.friendly_name, e.description
ORDER BY count DESC
