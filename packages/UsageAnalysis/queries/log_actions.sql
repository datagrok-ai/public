--name: Log Actions
--connection: System:Datagrok

select u.first_name, u.email, e.event_time, et.name, e.description
from events e
join event_types et on et.id = e.event_type_id
join users_sessions us on us.id = e.session_id
join users u on u.id = us.user_id
where u.first_name != 'System'
order by event_time desc
--end

--name: Log Actions Summary
--input: string eventTime = "today" {pattern: datetime}
--connection: System:Datagrok

select u.login, t.name, t.source, count(*) as runs from events e
inner join event_types t on t.id = e.event_type_id
inner join users_sessions s on s.id = e.session_id
inner join users u on u.id = s.user_id
where @eventTime(e.event_time)
group by u.login, t.name, t.source
order by 1,2,4 desc
--end

--name: Log Actions Summary by Hours
--input: string eventTime = "today" {pattern: datetime}
--connection: System:Datagrok

select u.login, t.name, t.source, date_part('hour', e.event_time) as hour, count(*) as runs from events e
inner join event_types t on t.id = e.event_type_id
inner join users_sessions s on s.id = e.session_id
inner join users u on u.id = s.user_id
where @eventTime(e.event_time)
group by u.login, t.name, t.source, hour
order by 1,2,4 desc
--end
