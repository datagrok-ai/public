--name: TopEventsOfUser
--meta.cache: true
--input: string date { pattern: datetime }
--input: string name
--connection: System:Datagrok
select et.friendly_name, count(1) from events e
inner join event_types et on e.event_type_id = et.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where
@date(e.event_time)
and u.login = @name
and et.friendly_name is not null
and et.friendly_name != ''
group by et.friendly_name
limit 50
--end

