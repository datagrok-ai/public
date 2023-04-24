--name:UniqueUsersOverview
--input: string date {pattern: datetime}
--input: list groups
--meta.invalidate: 0 * * * *
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
)
select date(e.event_time) as date, count(distinct u.id)
	from events e
	inner join users_sessions s on e.session_id = s.id
	inner join users u on u.id = s.user_id
	inner join selected_groups sg on u.group_id = sg.id
	WHERE @date(e.event_time)
group by date(e.event_time)
--end

--name:UniqueUsersStats
--input: string date {pattern: datetime}
--input: list groups
--meta.invalidate: 0 * * * *
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
)
select u.name as name, u.id as uid, count(e.id) as count
	from events e
	inner join users_sessions s on e.session_id = s.id
	inner join users u on u.id = s.user_id
	inner join selected_groups sg on u.group_id = sg.id
	WHERE @date(e.event_time)
group by u.name, u.id
--end


--name: PackagesUsageOverview
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
select e.id, e.event_time, u.friendly_name as user, u.id as uid, u.group_id as ugid,
coalesce(pp.name, pp1.name, pp2.name) as package, coalesce(pp.package_id, pp1.package_id, pp2.package_id) as pid
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

)
select coalesce(res.package, 'Core') as package, res.user, res.user as name, count(*) AS count,
res.uid, res.ugid, res.pid, max(event_time) as time_end, min(event_time) as time_start
from res, selected_groups sg
where res.ugid = sg.id
and (res.package = any(@packages) or @packages = ARRAY['all'])
GROUP BY res.package, res.user, res.uid, res.ugid, res.pid
--end
