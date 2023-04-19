--name: PackagesUsage
--input: string date {pattern: datetime}
--input: list groups
--input: list packages
--meta.cache: true
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
),
res AS (
select pp.name as package, u.friendly_name as user,
e.event_time as time_old, u.id as uid, u.group_id as ugid, pp.package_id as pid
from events e
inner join event_types et on e.event_type_id = et.id
inner join entities en on et.id = en.id
inner join published_packages pp on en.package_id = pp.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time)
union all
select pp.name as package, u.friendly_name as user,
e.event_time as time_old, u.id as uid, u.group_id as ugid, pp.package_id as pid
from events e
inner join event_parameter_values epv inner join event_parameters ep
on epv.parameter_id = ep.id and ep.name = 'package' on epv.event_id = e.id
inner join published_packages pp on pp.id::text = epv.value
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time)),
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
select res.package, res.user, count(*) AS count,
to_timestamp(floor((extract('epoch' from res.time_old) / trunc )) * trunc)
AT TIME ZONE 'UTC' as time_start,
to_timestamp(floor((extract('epoch' from res.time_old) / trunc )) * trunc)
AT TIME ZONE 'UTC' + trunc * interval '1 sec' as time_end,
res.uid, res.ugid, res.pid
from res, t2, selected_groups sg
where res.ugid = sg.id
and (res.package = any(@packages) or @packages = ARRAY['all'])
GROUP BY res.package, res.user, time_start, time_end, res.uid, res.ugid, res.pid
--end


--name: PackagesContextPaneFunctions
--input: int time_start
--input: int time_end
--input: string users
--input: string packages
--meta.cache: true
--connection: System:Datagrok
select pp.name as package, en.id, et.name, count(*)
from events e
inner join event_types et on e.event_type_id = et.id
inner join entities en on et.id = en.id
inner join published_packages pp on en.package_id = pp.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where e.event_time between to_timestamp(@time_start)
and to_timestamp(@time_end)
and u.id = any(@users)
and pp.package_id = any(@packages)
group by en.id, et.name, pp.name
--end


--name: PackagesContextPaneLogs
--input: int time_start
--input: int time_end
--input: string users
--input: string packages
--meta.cache: true
--connection: System:Datagrok
select et.source, count(*)
from events e
inner join event_types et on e.event_type_id = et.id
inner join event_parameter_values epv inner join event_parameters ep
on epv.parameter_id = ep.id and ep.name = 'package' and ep.type != 'entity_id'
on epv.event_id = e.id
inner join published_packages pp on pp.id::text = epv.value
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where e.event_time between to_timestamp(@time_start)
and to_timestamp(@time_end)
and u.id = any(@users)
and pp.package_id = any(@packages)
group by et.source
--end


--name: PackagesContextPaneAudit
--input: int time_start
--input: int time_end
--input: string users
--input: string packages
--meta.cache: true
--connection: System:Datagrok
select et.friendly_name as name, count(*)
from event_parameter_values e
inner join events ev on e.event_id = ev.id
inner join event_types et on ev.event_type_id = et.id
inner join event_parameters ep on ep.id = e.parameter_id and ep.type = 'entity_id'
inner join entities en on e.value != 'null' and en.id = e.value::uuid
left join published_packages pp inner join packages p1 on p1.id = pp.package_id on en.package_id = pp.id
left join packages p on p.id = en.id
inner join users_sessions s on ev.session_id = s.id
inner join users u on u.id = s.user_id
where (not p.id is null or not p1.id is null)
and ev.event_time between to_timestamp(@time_start)
and to_timestamp(@time_end)
and u.id = any(@users)
and COALESCE(p.id, p1.id) = any(@packages)
group by et.friendly_name
--end
