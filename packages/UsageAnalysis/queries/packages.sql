--name: PackagesUsage
--input: string date {pattern: datetime}
--input: list groups
--input: list packages
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
),
res AS (
select e.event_time as time_old, u.friendly_name as user, u.id as uid, u.group_id as ugid,
coalesce(pp.name, pp1.name, pp2.name, 'Core') as package,
coalesce(pp.package_id, pp1.package_id, pp2.package_id) as pid
from events e
inner join event_types et on e.event_type_id = et.id
left join entities en on et.id = en.id
left join published_packages pp on en.package_id = pp.id
left join event_parameter_values epv inner join event_parameters ep on epv.parameter_id = ep.id and ep.name = 'package'
on epv.event_id = e.id
left join published_packages pp1 on pp1.id::text = epv.value
left join event_parameter_values epv1 inner join event_parameters ep1 on epv1.parameter_id = ep1.id and ep1.type = 'entity_id'
inner join entities e1 on epv1.value != 'null' and e1.id = epv1.value::uuid
inner join published_packages pp2 inner join packages p2 on p2.id = pp2.package_id on e1.package_id = pp2.id
on epv1.event_id = e.id
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
select res.package, res.user, count(*) AS count,
to_timestamp(floor((extract('epoch' from res.time_old) / trunc )) * trunc)
AT TIME ZONE 'UTC' as time_start,
to_timestamp(floor((extract('epoch' from res.time_old) / trunc )) * trunc)
AT TIME ZONE 'UTC' + trunc * interval '1 sec' as time_end,
res.uid, res.ugid, coalesce(res.pid, '00000000-0000-0000-0000-000000000000') as pid
from res, t2, selected_groups sg
where res.ugid = sg.id
and (res.package = any(@packages) or @packages = ARRAY['all'])
GROUP BY res.package, res.user, time_start, time_end, res.uid, res.ugid, res.pid
--end


--name: PackagesContextPaneFunctions
--input: int time_start
--input: int time_end
--input: list users
--input: list packages
--connection: System:Datagrok
with res AS (
select coalesce(pp.name, p.name, 'Core') as package, en.id, et.name,
coalesce(pp.package_id, p1.id, '00000000-0000-0000-0000-000000000000') as pid
from events e
inner join event_types et on e.event_type_id = et.id
inner join entities en on et.id = en.id
left join project_relations pr ON pr.entity_id = en.id
left join projects p ON p.id = pr.project_id
and p.is_root = true
and p.is_package = true
left join published_packages pp on en.package_id = pp.id
left join packages p1 on p.name = p1.name or p.name = p1.friendly_name
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where e.event_time between to_timestamp(@time_start)
and to_timestamp(@time_end)
and u.id = any(@users)
and et.name is not null
)
select res.*, count(*)
from res
where res.pid = any(@packages)
group by res.package, res.id, res.name, res.pid
--end


--name: PackagesContextPaneLogs
--input: int time_start
--input: int time_end
--input: list users
--input: list packages
--connection: System:Datagrok
select et.source, count(*)
from events e
inner join event_types et on e.event_type_id = et.id
left join event_parameter_values epv inner join event_parameters ep
on epv.parameter_id = ep.id and ep.name = 'package' and ep.type != 'entity_id'
on epv.event_id = e.id
left join published_packages pp on pp.id::text = epv.value
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where et.source in ('debug', 'error', 'info', 'warning')
and e.event_time between to_timestamp(@time_start)
and to_timestamp(@time_end)
and u.id = any(@users)
and coalesce(pp.package_id, '00000000-0000-0000-0000-000000000000') = any(@packages)
group by et.source
--end


--name: PackagesContextPaneAudit
--input: int time_start
--input: int time_end
--input: list users
--input: list packages
--connection: System:Datagrok
select et.friendly_name as name, count(*)
from events e
inner join event_types et on e.event_type_id = et.id
inner join event_parameter_values epv on epv.event_id = e.id
inner join event_parameters ep on ep.id = epv.parameter_id and ep.type = 'entity_id'
inner join entities en on epv.value != 'null' and en.id::text = epv.value
left join published_packages pp inner join packages p1 on p1.id = pp.package_id on en.package_id = pp.id
left join packages p on p.id = en.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where et.source = 'audit'
and e.event_time between to_timestamp(@time_start)
and to_timestamp(@time_end)
and u.id = any(@users)
and coalesce(pp.package_id, p.id, '00000000-0000-0000-0000-000000000000') = any(@packages)
group by et.friendly_name
--end
