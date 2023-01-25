--name: NewUsersSummaryByDate
--input: string days
--connection: System:Datagrok
--meta.cache: true
--meta.invalidate: 0 0 * ? * * *
select t.date::date, count(t.id) as user_count from (
	select distinct on (date(e.event_time), u.id) u.id, e.event_time as date
	from events e
	inner join users_sessions s on e.session_id = s.id
	inner join users u on u.id = s.user_id
	WHERE DATE(e.event_time) BETWEEN (current_timestamp - (CONCAT(@days, ' day'))::interval)
    AND (current_timestamp - interval '0 day')
) t
group by t.date::date
order by t.date::date;
--end