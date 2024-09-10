--name: FunctionsUsage
--input: string date {pattern: datetime}
--input: list<string> groups
--input: list<string> packages
--meta.cache: all
--meta.cache.invalidateOn: 0 0 * * *
--connection: System:Datagrok
--test: FunctionsUsage(date='today', ['1ab8b38d-9c4e-4b1e-81c3-ae2bde3e12c5'], ['all'])
with recursive selected_groups as (
  select id from groups
  where id::varchar = any(@groups)
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
and (res.package = any(@packages) or @packages = ARRAY['all']::varchar[])
GROUP BY res.function, res.package, res.user, time_start, time_end,
res.uid, res.ugid, res.pid
--end


--name: FunctionsContextPane
--input: int time_start
--input: int time_end
--input: list<string> users
--input: list<string> packages
--input: list<string> functions
--meta.cache: all
--meta.cache.invalidateOn: 0 0 * * *
--connection: System:Datagrok
--test: FunctionsContextPane(1681084800, 1681516800, ['878c42b0-9a50-11e6-c537-6bf8e9ab02ee'], ['00000000-0000-0000-0000-000000000000'], ['OpenServerFile'])
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
and u.id::varchar = any(@users)
and et.name = any(@functions)
)
select res.package, res.run, res.function, res.time, res.rid, res.pid
from res
where res.pid::varchar = any(@packages)
--end


--name: FunctionsExecTime
--input: string function
--connection: System:Datagrok
with res as (
select e.id, EXTRACT(EPOCH FROM (e.event_time_finished - e.event_time)) AS time, ep.name as input,
ep.type as type, COALESCE(epv.value, epv.value_array::text, epv.data_frame_value_name, epv.value_string, epv.value_uuid::text) as value
from events e
inner join event_types et on e.event_type_id = et.id
left join event_parameter_values epv inner join event_parameters ep on epv.parameter_id = ep.id
and ep.is_input = true on epv.event_id = e.id
where et.name = @function)
select res.time, json_object_agg(res.input, res.value)::text as input
from res
GROUP BY res.time
--end