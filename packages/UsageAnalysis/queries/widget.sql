--name:UniqueUsersSummary
--meta.cache: all
--meta.cache.invalidateOn: 0 0 * * *
--connection: System:Datagrok
select date(e.event_time) as date, count(distinct u.id)
	from events e
	inner join users_sessions s on e.session_id = s.id
	inner join users u on u.id = s.user_id
	WHERE e.event_time > (current_timestamp - '60 day'::interval)
group by date(e.event_time)
--end

--name:UsersEventsSummary
--meta.cache: all
--meta.cache.invalidateOn: 0 0 * * *
--connection: System:Datagrok
select date(e.event_time) as date,  count(e.id)
	from events e
	inner join users_sessions s on e.session_id = s.id
	inner join users u on u.id = s.user_id
	WHERE e.event_time > (current_timestamp - '60 day'::interval)
group by date(e.event_time)
--end

--name:UsersErrorsSummary
--meta.cache: all
--meta.cache.invalidateOn: 0 0 * * *
--connection: System:Datagrok
select date(e.event_time) as date,  count(e.id)
	from events e
	inner join event_types t on t.id = e.event_type_id
	inner join users_sessions s on e.session_id = s.id
	inner join users u on u.id = s.user_id
	WHERE e.event_time > (current_timestamp - '60 day'::interval) and t.source = 'error'
group by date(e.event_time)
--end
