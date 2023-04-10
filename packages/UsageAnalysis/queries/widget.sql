--name: WidgetData
--meta.cache: true
--connection: System:Datagrok
select COALESCE(u.name, 'System') as name,
to_timestamp(floor((extract('epoch' from e.event_time) / 3600 )) * 3600)
AT TIME ZONE 'UTC' as time, count(*)
from events e
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where e.event_time::TIMESTAMP between 'now'::TIMESTAMP - interval '7 day' and 'now'::TIMESTAMP
group by time, u.name
--end
