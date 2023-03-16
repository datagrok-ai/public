--name: PackagesUsage
--input: string date {pattern: datetime}
--input: list users
--input: list packages
--meta.cache: true
--connection: System:Datagrok
select t.package, t.user, sum(t.count)::bigint as count, t.time, t.login
from (
(select pp.name as package, u.friendly_name as user, count(*),
to_timestamp(floor((extract('epoch' from e.event_time) / 600 )) * 600)
AT TIME ZONE 'UTC' as time, u.login as login
from events e
inner join event_types et on e.event_type_id = et.id
inner join entities en on et.id = en.id
inner join published_packages pp on en.package_id = pp.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time)
and (u.login = any(@users) or @users = ARRAY['all'])
and (pp.name = any(@packages) or @packages = ARRAY['all'])
group by pp.name, u.friendly_name, time, u.login
limit 100000)
union all
(select pp.name as package, u.friendly_name as user, count(*),
to_timestamp(floor((extract('epoch' from e.event_time) / 600 )) * 600)
AT TIME ZONE 'UTC' as time, u.login as login
from events e
inner join event_parameter_values epv inner join event_parameters ep on epv.parameter_id = ep.id and ep.name = 'package' on epv.event_id = e.id
inner join published_packages pp on pp.id::text = epv.value
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time)
and (u.login = any(@users) or @users = ARRAY['all'])
and (pp.name = any(@packages) or @packages = ARRAY['all'])
group by pp.name, u.friendly_name, time, u.login
limit 100000)
) t
group by t.package, t.user, t.time, t.login
--end


--name: PackageUsageAtPoint
--input: int date
--input: string user
--input: string package
--meta.cache: true
--connection: System:Datagrok
select t.package, t.user, t.source, t.name, t.time
from (
(select pp.name as package, u.friendly_name as user, et.source, et.friendly_name as name, e.event_time as time
from events e
inner join event_types et on e.event_type_id = et.id
inner join entities en on et.id = en.id
inner join published_packages pp on en.package_id = pp.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where e.event_time::TIMESTAMP between to_timestamp(@date) and to_timestamp(@date + 600)
and u.login = @user
and pp.name = @package
limit 100000)
union all
(select pp.name as package, u.friendly_name as user, et.source, et.friendly_name as name, e.event_time as time
from events e
inner join event_types et on e.event_type_id = et.id
inner join event_parameter_values epv inner join event_parameters ep on epv.parameter_id = ep.id and ep.name = 'package' on epv.event_id = e.id
inner join published_packages pp on pp.id::text = epv.value
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where e.event_time::TIMESTAMP between to_timestamp(@date) and to_timestamp(@date + 600)
and u.login = @user
and pp.name = @package
limit 100000)
) t
--end
