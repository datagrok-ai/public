--name:UniqueUsersSummary
--meta.cache: all
--meta.cache.invalidateOn: 0 0 * * *
--connection: System:Datagrok
select date, count(*) from (
	select distinct date(e.event_time) as date, s.user_id
	from events e
	inner join users_sessions s on e.session_id = s.id
	WHERE e.event_time > (current_timestamp - '60 day'::interval)
) sub
group by date
--end

--name:UsersEventsSummary
--meta.cache: all
--meta.cache.invalidateOn: 0 0 * * *
--connection: System:Datagrok
select date(e.event_time) as date, count(e.id)
	from events e
	WHERE e.event_time > (current_timestamp - '60 day'::interval)
	  AND e.session_id IS NOT NULL
group by date(e.event_time)
--end

--name:UsersErrorsSummary
--meta.cache: all
--meta.cache.invalidateOn: 0 0 * * *
--connection: System:Datagrok
select date(e.event_time) as date, count(e.id)
	from events e
	WHERE e.event_time > (current_timestamp - '60 day'::interval)
	  AND e.session_id IS NOT NULL
	  AND e.event_type_id IN (SELECT id FROM event_types WHERE source = 'error')
group by date(e.event_time)
--end
