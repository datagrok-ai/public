--name: MostRecentEntities
--connection: System:Datagrok
--input: string user
WITH user_events AS (
  SELECT e.id, e.event_time
  FROM events e
  JOIN event_parameter_values vu ON vu.event_id = e.id
  JOIN event_parameters pu ON pu.id = vu.parameter_id
  WHERE pu.name = 'user'
    AND vu.value_uuid = @user
)
SELECT
  en.id,
  MAX(ue.event_time) AS last_event_time
FROM user_events ue
JOIN event_parameter_values v ON v.event_id = ue.id
JOIN event_parameters p ON p.id = v.parameter_id
JOIN entities en ON en.id = v.value_uuid
WHERE p.type = 'entity_id'
GROUP BY en.id
ORDER BY last_event_time DESC
LIMIT 50;
--end
