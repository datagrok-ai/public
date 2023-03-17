--name: PackagesUsage
--input: string date {pattern: datetime}
--input: list users
--input: list packages
--meta.cache: true
--connection: System:Datagrok
WITH res AS (
(select pp.name as package, u.friendly_name as user,
e.event_time as time_old, u.id as uid, u.group_id as ugid
from events e
inner join event_types et on e.event_type_id = et.id
inner join entities en on et.id = en.id
inner join published_packages pp on en.package_id = pp.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time)
and (u.login = any(@users) or @users = ARRAY['all'])
and (pp.name = any(@packages) or @packages = ARRAY['all'])
limit 100000)
union all
(select pp.name as package, u.friendly_name as user,
e.event_time as time_old, u.id as uid, u.group_id as ugid
from events e
inner join event_parameter_values epv inner join event_parameters ep
on epv.parameter_id = ep.id and ep.name = 'package' on epv.event_id = e.id
inner join published_packages pp on pp.id::text = epv.value
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time)
and (u.login = any(@users) or @users = ARRAY['all'])
and (pp.name = any(@packages) or @packages = ARRAY['all'])
limit 100000)
),
t1 AS (
  SELECT (MAX(res.time_old) - MIN(res.time_old)) as inter
  FROM res
),
t2 AS (
  SELECT case when inter >= INTERVAL '6 month' then 432000
    when inter >= INTERVAL '2 month' then 144000
	when inter >= INTERVAL '1 week' then 86400
	when inter >= interval '1 day' then 14400
	else 600 end as trunc
from t1
)
select res.package, res.user, count(*) AS count,
to_timestamp(floor((extract('epoch' from res.time_old) / trunc )) * trunc)
AT TIME ZONE 'UTC' as time, res.uid, res.ugid
from res, t2
GROUP BY res.package, res.user, time, res.uid, res.ugid
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
