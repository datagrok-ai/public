--name:UniqueUsersSummary
--meta.cache: true
--meta.invalidate: 0 * * * *
--connection: System:Datagrok
select date(e.event_time) as date, count(distinct u.id)
	from events e
	inner join users_sessions s on e.session_id = s.id
	inner join users u on u.id = s.user_id
	WHERE e.event_time > (current_timestamp - '60 day'::interval)
group by date(e.event_time)
--end

--name:UsersEventsSummary
--meta.cache: true
--meta.invalidate: 0 * * * *
--connection: System:Datagrok
select date(e.event_time) as date,  count(e.id)
	from events e
	inner join users_sessions s on e.session_id = s.id
	inner join users u on u.id = s.user_id
	WHERE e.event_time > (current_timestamp - '60 day'::interval)
group by date(e.event_time)
--end
