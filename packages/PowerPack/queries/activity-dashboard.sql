--name: MostRecentEntities
--connection: System:Datagrok
--input: string user
WITH user_events AS (
  SELECT e.id, e.event_time
  FROM events e
  JOIN event_parameter_values vu ON vu.event_id = e.id
  JOIN event_parameters pu ON pu.id = vu.parameter_id
  WHERE pu.name = 'user' AND vu.value_uuid = @user
    AND e.event_time > NOW() - INTERVAL '30 days'
)
SELECT
  en.id,
  EXTRACT(EPOCH FROM MAX(ue.event_time) AT TIME ZONE 'UTC') * 1000 AS last_event_time
FROM user_events ue
JOIN event_parameter_values v ON v.event_id = ue.id
JOIN event_parameters p ON p.id = v.parameter_id
JOIN entities en ON en.id = v.value_uuid
JOIN entity_types et ON et.id = en.entity_type_id
WHERE p.type = 'entity_id'
  AND et.name NOT IN ('FuncCall', 'UserGroup', 'User', 'GrokPackage', 'UserReport',
                      'TableInfo', 'PackageFunc', 'Func', 'ViewInfo', 'DataConnection',
                      'ColumnInfo')
GROUP BY en.id
ORDER BY last_event_time DESC
LIMIT 16;
--end

--name: RecentlySharedWithMe
--connection: System:Datagrok
--input: string user
WITH RECURSIVE user_groups AS (
  SELECT group_id AS gid FROM users WHERE id = @user
  UNION
  SELECT r.parent_id FROM groups_relations r
  JOIN user_groups ug ON r.child_id = ug.gid
),
share_events AS (
  SELECT e.id, e.event_time
  FROM events e
  JOIN event_types et ON et.id = e.event_type_id
  JOIN event_parameter_values vu ON vu.event_id = e.id
  JOIN event_parameters pu ON pu.id = vu.parameter_id
  WHERE et.friendly_name IN ('entity-shared', 'entity-shared-silent')
    AND (
      (pu.name = 'grantedUser' AND vu.value_uuid = @user)
      OR (pu.name = 'group' AND vu.value_uuid IN (SELECT gid FROM user_groups))
    )
    AND e.event_time > NOW() - INTERVAL '30 days'
)
SELECT
  en.id,
  EXTRACT(EPOCH FROM MAX(se.event_time) AT TIME ZONE 'UTC') * 1000 AS shared_time
FROM share_events se
JOIN event_parameter_values v ON v.event_id = se.id
JOIN event_parameters p ON p.id = v.parameter_id
JOIN entities en ON en.id = v.value_uuid
JOIN entity_types et ON et.id = en.entity_type_id
WHERE p.name = 'entity'
  AND et.name NOT IN ('FuncCall', 'UserGroup', 'User', 'GrokPackage', 'UserReport',
                      'TableInfo', 'PackageFunc', 'Func', 'ViewInfo', 'DataConnection',
                      'ColumnInfo')
GROUP BY en.id
ORDER BY shared_time DESC
LIMIT 8;
--end
