--name: FunctionsUsage
--input: string date {pattern: datetime}
--input: list groups
--input: list packages
--meta.cache: true
--meta.invalidate: 0 * * * *
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
),
res AS (
select DISTINCT e.id as id_, et.name as function, coalesce(pp.name, p1.name, 'Core') as package,
u.friendly_name as user, e.event_time as time_old, u.id as uid,
u.group_id as ugid, coalesce(pp.package_id, p1.id, '00000000-0000-0000-0000-000000000000') as pid
from events e
inner join event_types et on e.event_type_id = et.id
inner join entities en on et.id = en.id
left join published_packages pp on en.package_id = pp.id
left join project_relations pr ON pr.entity_id = en.id
left join projects proj ON proj.id = pr.project_id
and proj.is_root = true
and proj.is_package = true
left join packages p1 on proj.name = p1.name or proj.name = p1.friendly_name
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time)
),
t1 AS (
  SELECT (MAX(res.time_old) - MIN(res.time_old)) as inter
  FROM res
),
t2 AS (
  SELECT case when inter >= INTERVAL '6 month' then 864000
  when inter >= INTERVAL '70 day' then 432000
  when inter >= INTERVAL '10 day' then 86400
  when inter >= INTERVAL '2 day' then 21600
  when inter >= INTERVAL '3 hour' then 3600
  else 600 end as trunc
from t1
)
select res.function, res.package, res.user, count(*) AS count,
to_timestamp(floor((extract('epoch' from res.time_old) / trunc )) * trunc)
AT TIME ZONE 'UTC' as time_start,
to_timestamp(floor((extract('epoch' from res.time_old) / trunc )) * trunc)
AT TIME ZONE 'UTC' + trunc * interval '1 sec' as time_end,
res.uid, res.ugid, coalesce(res.pid, '00000000-0000-0000-0000-000000000000') as pid
from res, t2, selected_groups sg
where res.ugid = sg.id
and (res.package = any(@packages) or @packages = ARRAY['all'])
GROUP BY res.function, res.package, res.user, time_start, time_end,
res.uid, res.ugid, res.pid
--end


--name: FunctionsContextPane
--input: int time_start
--input: int time_end
--input: list users
--input: list packages
--input: list functions
--meta.cache: true
--meta.invalidate: 0 * * * *
--connection: System:Datagrok
with res AS (
select DISTINCT e.id as id_, coalesce(pp.name, p1.name, 'Core') as package,
e.friendly_name as run, et.name as function, e.event_time as time, e.id as rid,
coalesce(pp.package_id, p1.id, '00000000-0000-0000-0000-000000000000') as pid
from events e
inner join event_types et on e.event_type_id = et.id
inner join entities en on et.id = en.id
left join published_packages pp on en.package_id = pp.id
left join project_relations pr ON pr.entity_id = en.id
left join projects proj ON proj.id = pr.project_id
and proj.is_root = true
and proj.is_package = true
left join packages p1 on proj.name = p1.name or proj.name = p1.friendly_name
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where e.event_time between to_timestamp(@time_start)
and to_timestamp(@time_end)
and u.id = any(@users)
and et.name = any(@functions)
)
select res.package, res.run, res.function, res.time, res.rid, res.pid
from res
where res.pid = any(@packages)
--end
