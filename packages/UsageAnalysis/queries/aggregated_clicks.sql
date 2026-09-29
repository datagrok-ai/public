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
  AND et.friendly_name IN ('click', 'menu click', 'dialog show', 'dialog close', 'dialog ok', 'input', 'command', 'navigate', 'dragstart', 'drop', 'dock', 'undock')
  AND @date(e.event_time)
GROUP BY et.friendly_name, e.description
ORDER BY count DESC
--end


--name: Clicks
--friendlyName: Clicks
--input: string date {pattern: datetime}
--input: list<string> groups
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
)
select e.event_time, u.friendly_name as user, et.friendly_name as event_type, e.description, e.request_id,
  u.group_id as ugid, e.id
from events e
join event_types et on e.event_type_id = et.id
left join users_sessions s on e.session_id = s.id
left join users u on u.id = s.user_id
where et.source = 'usage'
  and et.friendly_name in ('click', 'menu click', 'dialog show', 'dialog close', 'dialog ok', 'input', 'command', 'navigate', 'dragstart', 'drop', 'dock', 'undock')
  and @date(e.event_time)
  and (u.id is null or u.group_id in (select id from selected_groups))
order by e.event_time desc
limit 10000
--end
