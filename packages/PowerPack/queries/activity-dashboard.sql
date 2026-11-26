--name: MostRecentEntities
--friendlyName: Most Recent Entities
--description: Entities the user interacted with most recently (last 30 days)
--connection: System:Datagrok
--input: string user [User ID whose recent activity to fetch]
-- Looks at the latest 5000 events of the user and never reads all of them: an account with
-- hundreds of thousands of events a month (a test user) otherwise loaded the database on every call.
-- The query picks the way itself, so it does not depend on the plan Postgres caches for it:
-- up to 20000 values with the user's id (n): start from those values;
-- more: walk the 200000 newest events of the instance and keep the user's.
-- OFFSET 0 keeps each subquery as written.
WITH n AS MATERIALIZED (
  SELECT count(*) AS cnt
  FROM (SELECT 1 FROM event_parameter_values WHERE value_uuid = @user LIMIT 20001) x
), user_events AS MATERIALIZED (
  (SELECT e.id, e.event_type_id, e.event_time
   FROM (
     SELECT vu.event_id
     FROM event_parameter_values vu
     JOIN event_parameters pu ON pu.id = vu.parameter_id AND pu.name = 'user'
     WHERE vu.value_uuid = @user AND (SELECT cnt FROM n) <= 20000
     OFFSET 0) u
   JOIN events e ON e.id = u.event_id
   WHERE e.event_time > NOW() - INTERVAL '30 days'
   ORDER BY e.event_time DESC
   LIMIT 5000)
  UNION ALL
  (SELECT r.id, r.event_type_id, r.event_time
   FROM (
     SELECT e.id, e.event_type_id, e.event_time
     FROM events e
     WHERE e.event_time > NOW() - INTERVAL '30 days'
     ORDER BY e.event_time DESC
     LIMIT 200000) r
   WHERE (SELECT cnt FROM n) > 20000
     AND EXISTS (
       SELECT 1
       FROM event_parameter_values vu
       JOIN event_parameters pu ON pu.id = vu.parameter_id AND pu.name = 'user'
       WHERE vu.event_id = r.id AND vu.value_uuid = @user
       OFFSET 0)
   ORDER BY r.event_time DESC
   LIMIT 5000)
)
SELECT
  en.id,
  EXTRACT(EPOCH FROM MAX(ue.event_time) AT TIME ZONE 'UTC') * 1000 AS last_event_time
FROM user_events ue
JOIN event_parameters p ON p.event_type_id = ue.event_type_id AND p.type = 'entity_id'
JOIN event_parameter_values v ON v.event_id = ue.id AND v.parameter_id = p.id
JOIN entities en ON en.id = v.value_uuid
JOIN entity_types et ON et.id = en.entity_type_id
WHERE et.name <> ALL (ARRAY['FuncCall', 'UserGroup', 'User', 'GrokPackage', 'UserReport',
    'TableInfo', 'PackageFunc', 'Func', 'ViewInfo', 'DataConnection', 'ColumnInfo'])
GROUP BY en.id
ORDER BY last_event_time DESC
LIMIT 16;
--end

--name: RecentlySharedWithMe
--friendlyName: Recently Shared With Me
--description: Entities recently shared with the user or their groups (last 30 days)
--connection: System:Datagrok
--input: string user [User ID whose shared entities to fetch]
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
