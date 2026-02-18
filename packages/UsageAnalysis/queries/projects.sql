--name: UniqueUsersPerProject
--input: string date {pattern: datetime}
--input: list<string> groups
--input: list<string> projects
--meta.cache: all
--meta.cache.invalidateOn: 0 0 * * *
--connection: System:Datagrok
WITH RECURSIVE selected_groups AS
(
    SELECT id FROM groups
    WHERE id = ANY(@groups)
    UNION
    SELECT gr.child_id AS id FROM selected_groups sg
    JOIN groups_relations gr ON sg.id = gr.parent_id
)
SELECT COUNT(DISTINCT u.id) AS unique_users, COALESCE(NULLIF(p.name, ''), p.friendly_name) AS project_name
FROM events e
INNER JOIN event_types et ON et.id = e.event_type_id AND et.source = 'usage' AND et.friendly_name = 'project-opened'
INNER JOIN users_sessions s ON e.session_id = s.id
INNER JOIN users u ON u.id = s.user_id
INNER JOIN event_parameter_values evv ON evv.value_uuid IS NOT NULL AND evv.event_id = e.id
INNER JOIN projects p ON p.id = evv.value_uuid
JOIN selected_groups sg ON u.group_id = sg.id
WHERE @date(e.event_time)
AND p.is_dashboard AND NOT p.is_root AND NOT p.is_entity
AND COALESCE(p.name, p.friendly_name, '') <> ''
AND (@projects = ARRAY['all'] OR COALESCE(NULLIF(p.name, ''), p.friendly_name) = ANY(@projects))
GROUP BY COALESCE(NULLIF(p.name, ''), p.friendly_name)
ORDER BY unique_users DESC
--end

--name: DailyProjectAccess
--input: string date {pattern: datetime}
--input: list<string> groups
--input: list<string> projects
--meta.cache: all
--meta.cache.invalidateOn: 0 0 * * *
--connection: System:Datagrok
WITH RECURSIVE selected_groups AS
(
    SELECT id FROM groups
    WHERE id = ANY(@groups)
    UNION
    SELECT gr.child_id AS id FROM selected_groups sg
    JOIN groups_relations gr ON sg.id = gr.parent_id
)
SELECT
    COALESCE(NULLIF(p.name, ''), p.friendly_name) AS project_name,
    DATE(e.event_time) AS access_date,
    u.friendly_name AS user_name
FROM events e
INNER JOIN event_types et ON et.id = e.event_type_id AND et.source = 'usage' AND et.friendly_name = 'project-opened'
INNER JOIN users_sessions s ON e.session_id = s.id
INNER JOIN users u ON u.id = s.user_id
INNER JOIN event_parameter_values evv ON evv.value_uuid IS NOT NULL AND evv.event_id = e.id
INNER JOIN projects p ON p.id = evv.value_uuid
JOIN selected_groups sg ON u.group_id = sg.id
WHERE @date(e.event_time)
AND p.is_dashboard AND NOT p.is_root AND NOT p.is_entity
AND COALESCE(p.name, p.friendly_name, '') <> ''
AND (@projects = ARRAY['all'] OR COALESCE(NULLIF(p.name, ''), p.friendly_name) = ANY(@projects))
GROUP BY COALESCE(NULLIF(p.name, ''), p.friendly_name), DATE(e.event_time), u.friendly_name
ORDER BY access_date DESC, project_name, user_name
--end

--name: AccessCountPerPeriodPerProject
--input: string date {pattern: datetime}
--input: list<string> groups
--input: list<string> projects
--meta.cache: all
--meta.cache.invalidateOn: 0 0 * * *
--connection: System:Datagrok
WITH RECURSIVE selected_groups AS
(
    SELECT id FROM groups
    WHERE id = ANY(@groups)
    UNION
    SELECT gr.child_id AS id FROM selected_groups sg
    JOIN groups_relations gr ON sg.id = gr.parent_id
),
filtered_events AS (
    SELECT e.id AS event_id, e.event_time, COALESCE(NULLIF(p.name, ''), p.friendly_name) AS project_name
    FROM events e
    INNER JOIN event_types et ON et.id = e.event_type_id AND et.source = 'usage' AND et.friendly_name = 'project-opened'
    INNER JOIN users_sessions s ON e.session_id = s.id
    INNER JOIN users u ON u.id = s.user_id
    INNER JOIN event_parameter_values evv ON evv.value_uuid IS NOT NULL AND evv.event_id = e.id
    INNER JOIN projects p ON p.id = evv.value_uuid
    JOIN selected_groups sg ON u.group_id = sg.id
    WHERE @date(e.event_time)
    AND p.is_dashboard AND NOT p.is_root AND NOT p.is_entity
    AND COALESCE(p.name, p.friendly_name, '') <> ''
    AND (@projects = ARRAY['all'] OR COALESCE(NULLIF(p.name, ''), p.friendly_name) = ANY(@projects))
),
date_range AS (
    SELECT
        MIN(event_time) AS min_date,
        MAX(event_time) AS max_date,
        MAX(event_time) - MIN(event_time) AS date_diff
    FROM filtered_events
),
aggregated_events AS (
    SELECT fe.project_name,
        CASE
        WHEN dr.date_diff <= INTERVAL '31 days'
        THEN DATE_TRUNC('day', fe.event_time)
        WHEN dr.date_diff <= INTERVAL '1 year'
        THEN DATE_TRUNC('month', fe.event_time)
        ELSE DATE_TRUNC('year', fe.event_time)
        END AS period,
        COUNT(*) AS access_count
    FROM filtered_events fe
    CROSS JOIN date_range dr
    GROUP BY fe.project_name, period
    ORDER BY fe.project_name, period
)
SELECT * FROM aggregated_events;
--end

--name: ProjectsList
--connection: System:Datagrok
SELECT DISTINCT COALESCE(NULLIF(p.name, ''), p.friendly_name) AS project_name
FROM projects p
WHERE p.is_dashboard AND NOT p.is_root AND NOT p.is_entity AND COALESCE (p.name, p.friendly_name, '') <> '';
--end
