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
),
res AS
(
    SELECT COUNT(DISTINCT u.id) AS unique_users, COALESCE(NULLIF(p.name, ''), p.friendly_name) AS project_name, pp.package_id AS package_id, u.group_id AS ugid FROM events e
    INNER JOIN event_types et ON et.id = e.event_type_id AND et.source = 'usage' AND et.friendly_name = 'project-opened'
    INNER JOIN users_sessions s ON e.session_id = s.id
    INNER JOIN users u ON u.id = s.user_id
    INNER JOIN event_parameter_values evv ON evv.value_uuid IS NOT NULL AND evv.event_id = e.id
    INNER JOIN projects p ON p.id = evv.value_uuid
    INNER JOIN entities en ON p.id = en.id
    LEFT JOIN published_packages pp ON en.package_id = pp.id
    WHERE @date(e.event_time)
    AND p.is_dashboard AND NOT p.is_root AND NOT p.is_entity
    AND COALESCE (p.name, p.friendly_name, '') <> ''
    GROUP BY p.name, p.friendly_name, pp.package_id, u.group_id
),
res1 as (
    SELECT res.unique_users, res.project_name, res.package_id FROM res, selected_groups
    WHERE res.ugid = selected_groups.id
    AND (@projects = ARRAY['all'] OR res.project_name = ANY(@projects))
    ORDER BY res.unique_users DESC
)
select (sum(res1.unique_users)::int) unique_users, res1.project_name from res1
group by res1.project_name, res1.package_id
--end

--name: UserAccessFrequencyPerProject
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
user_access AS (
    SELECT
        u.id AS user_id,
        u.group_id AS ugid,
        pp.package_id as package_id,
        COALESCE(NULLIF(p.name, ''), p.friendly_name) AS project_name,
        DATE(e.event_time) AS access_date
    FROM events e
    INNER JOIN event_types et ON et.id = e.event_type_id AND et.source = 'usage' AND et.friendly_name = 'project-opened'
    INNER JOIN users_sessions s ON e.session_id = s.id
    INNER JOIN users u ON u.id = s.user_id
    INNER JOIN event_parameter_values evv ON evv.value_uuid IS NOT NULL AND evv.event_id = e.id
    INNER JOIN projects p ON p.id = evv.value_uuid
    INNER JOIN entities en ON p.id = en.id
    LEFT JOIN published_packages pp ON en.package_id = pp.id
    WHERE @date(e.event_time)
    AND p.is_dashboard AND NOT p.is_root AND NOT p.is_entity
    AND COALESCE (p.name, p.friendly_name, '') <> ''
    AND (@projects = ARRAY['all'] OR COALESCE(NULLIF(p.name, ''), p.friendly_name) = ANY(@projects))
    GROUP BY u.id, u.group_id, pp.id, p.name, p.friendly_name, DATE(e.event_time)
)
SELECT
    ua.project_name,
    (SELECT COUNT(DISTINCT user_id) FROM user_access WHERE project_name = ua.project_name) AS unique_users,
    ROUND(AVG(access_gap), 2) AS avg_days_between_accesses
FROM (
         SELECT
             user_id,
             package_id,
             project_name,
             access_date,
             ugid,
             LEAD(access_date) OVER (PARTITION BY user_id, project_name, package_id ORDER BY access_date) - access_date AS access_gap
         FROM user_access
     ) subquery
         JOIN user_access ua ON subquery.project_name = ua.project_name
WHERE access_gap IS NOT NULL
GROUP BY ua.package_id, ua.project_name
ORDER BY avg_days_between_accesses, unique_users desc
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
   SELECT e.id AS event_id, e.event_time, COALESCE(NULLIF(p.name, ''), p.friendly_name) AS project_name, u.group_id AS ugid, pp.package_id
   FROM events e
INNER JOIN event_types et ON et.id = e.event_type_id AND et.source = 'usage' AND et.friendly_name = 'project-opened'
INNER JOIN users_sessions s ON e.session_id = s.id
INNER JOIN users u ON u.id = s.user_id
INNER JOIN event_parameter_values evv ON evv.value_uuid IS NOT NULL AND evv.event_id = e.id
INNER JOIN projects p ON p.id = evv.value_uuid
INNER JOIN entities en ON p.id = en.id
LEFT JOIN published_packages pp ON en.package_id = pp.id
   WHERE @date(e.event_time)
   AND p.is_dashboard AND NOT p.is_root AND NOT p.is_entity
   AND COALESCE (p.name, p.friendly_name, '') <> ''
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
SELECT fe.project_name, fe.package_id,
    CASE
    WHEN dr.date_diff <= INTERVAL '31 days'
    THEN DATE_TRUNC('day', fe.event_time)
    WHEN dr.date_diff <= INTERVAL '1 year'
    THEN DATE_TRUNC('month', fe.event_time)
    ELSE DATE_TRUNC('year', fe.event_time)
    END AS period,
    COUNT(*) AS access_count
FROM filtered_events fe
    JOIN selected_groups sg ON sg.id = fe.ugid
    CROSS JOIN date_range dr
GROUP BY fe.project_name, period, fe.package_id
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
