--name: FunctionsUsage
--input: string date {pattern: datetime}
--input: list users
--input: list packages
--meta.cache: true
--connection: System:Datagrok
select et.friendly_name as function, pp.name as package, u.friendly_name as user, count(*),
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
and (pp.name = any(@packages) or @packages = ARRAY['all'])
group by et.friendly_name, pp.name, u.friendly_name, time
limit 10000
--end
