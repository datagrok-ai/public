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
select e.event_time as time, u.friendly_name as user, et.friendly_name as type, e.description as element,
  e.request_id as "action id"
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

--name: ClicksFollowedByError
--friendlyName: Clicks Followed By Error
--input: string date {pattern: datetime}
--input: list<string> groups
--input: string view = ""
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
), clicks as (
  select e.event_time, e.description, e.request_id,
    coalesce(u.id::text, (select v.value from event_parameter_values v join event_parameters p on p.id = v.parameter_id
      where v.event_id = e.id and p.name = 'anonSession' limit 1)) as who
  from events e
  join event_types et on e.event_type_id = et.id
  left join users_sessions s on e.session_id = s.id
  left join users u on u.id = s.user_id
  where et.source = 'usage'
    and et.friendly_name in ('click', 'menu click')
    and @date(e.event_time)
    and (u.id is null or u.group_id in (select id from selected_groups))
    and (coalesce(@view, '') = '' or exists (select 1 from event_parameter_values v
      join event_parameters p on p.id = v.parameter_id
      where v.event_id = e.id and p.name = 'view' and lower(v.value) = lower(@view)))
), marked as (
  select c.*, exists (select 1 from events x join event_types xt on xt.id = x.event_type_id
    where xt.source = 'error' and x.event_time >= c.event_time and x.event_time <= c.event_time + interval '5 seconds'
      and (x.request_id = c.request_id or x.request_id like c.request_id || '.%')) as followed
  from clicks c
)
select m.description as element, count(*)::int as clicks, count(distinct m.who)::int as users,
  (count(*) filter (where m.followed))::int as followed_by_error,
  round(100.0 * count(*) filter (where m.followed) / count(*), 1)::float as followed_by_error_pct
from marked m
group by m.description
order by followed_by_error desc, clicks desc
limit 1000
--end

--name: ClickErrors
--friendlyName: Click Errors
--input: string date {pattern: datetime}
--input: list<string> groups
--input: string element
--input: string view = ""
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
), clicks as (
  select e.event_time, e.request_id
  from events e
  join event_types et on e.event_type_id = et.id
  left join users_sessions s on e.session_id = s.id
  left join users u on u.id = s.user_id
  where et.source = 'usage'
    and et.friendly_name in ('click', 'menu click')
    and e.description = @element
    and e.request_id is not null
    and @date(e.event_time)
    and (u.id is null or u.group_id in (select id from selected_groups))
    and (coalesce(@view, '') = '' or exists (select 1 from event_parameter_values v
      join event_parameters p on p.id = v.parameter_id
      where v.event_id = e.id and p.name = 'view' and lower(v.value) = lower(@view)))
)
select t.error_stack_trace_hash::text as signature, min(coalesce(nullif(x.error_message, ''), x.description)) as error,
  count(distinct c.request_id)::int as clicks, max(c.request_id) as action
from clicks c
join events x on x.event_time >= c.event_time and x.event_time <= c.event_time + interval '5 seconds'
  and (x.request_id = c.request_id or x.request_id like c.request_id || '.%')
join event_types t on t.id = x.event_type_id and t.source = 'error'
group by 1
order by clicks desc
limit 10
--end