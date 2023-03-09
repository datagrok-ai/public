--name: PackagesUsage
--input: string date {pattern: datetime}
--input: list users
--meta.cache: true
--connection: System:Datagrok
select t.package, t.user, sum(t.count)::bigint as count, t.time
from (
(select pp.name as package, u.friendly_name as user, count(*),
to_timestamp(floor((extract('epoch' from e.event_time) / 600 )) * 600)
AT TIME ZONE 'UTC' as time
from events e
inner join event_types et on e.event_type_id = et.id
inner join entities en on et.id = en.id
inner join published_packages pp on en.package_id = pp.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time)
and (u.login = any(@users) or @users = ARRAY['all'])
-- where e.event_time::TIMESTAMP between 'now'::TIMESTAMP - interval '7 day' and 'now'::TIMESTAMP
group by pp.name, u.friendly_name, time
limit 100000)
union all
(select pp.name as package, u.friendly_name as user, count(*),
to_timestamp(floor((extract('epoch' from e.event_time) / 600 )) * 600)
AT TIME ZONE 'UTC' as time
from events e
inner join event_parameter_values epv inner join event_parameters ep on epv.parameter_id = ep.id and ep.name = 'package' on epv.event_id = e.id
inner join published_packages pp on pp.id::text = epv.value
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time)
and (u.login = any(@users) or @users = ARRAY['all'])
-- where e.event_time::TIMESTAMP between 'now'::TIMESTAMP - interval '7 day' and 'now'::TIMESTAMP
group by pp.name, u.friendly_name, time
limit 100000)
) t
group by t.package, t.user, t.time
--end
