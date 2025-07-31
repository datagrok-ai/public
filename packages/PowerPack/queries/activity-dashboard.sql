--name: MostRecentEntities
--connection: System:Datagrok
--input: string user
SELECT
  en.id,
  MAX(e.event_time) AS last_event_time
FROM events e
JOIN event_parameter_values v ON v.event_id = e.id
JOIN event_parameters p ON p.id = v.parameter_id AND p.type = 'entity_id'
JOIN entities en ON en.id = v.value_uuid
JOIN event_parameter_values vu ON vu.event_id = e.id
JOIN event_parameters pu ON pu.id = vu.parameter_id AND pu.name = 'user'
WHERE
  vu.value_uuid = @user
GROUP BY en.id
ORDER BY last_event_time DESC
LIMIT 20
--end
