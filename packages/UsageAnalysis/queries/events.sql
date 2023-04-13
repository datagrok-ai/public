--name: EventsSources
--meta.cache: true
--input: string date {pattern: datetime}
--connection: System:Datagrok
with res as (
select et.source, e.event_time as time_old
from events e
inner join event_types et on e.event_type_id = et.id
where @date(e.event_time)),
t1 AS (
  SELECT (MAX(res.time_old) - MIN(res.time_old)) as inter
  FROM res
),
t2 AS (
  SELECT case when inter >= INTERVAL '6 month' then 864000
  when inter >= INTERVAL '70 day' then 432000
	when inter >= INTERVAL '10 day' then 86400
	when inter >= INTERVAL '2 day' then 21600
	when inter >= INTERVAL '3 hour' then 3600
	else 600 end as trunc
from t1
)
select res.source, count(*),
to_timestamp(floor((extract('epoch' from res.time_old) / trunc )) * trunc)
AT TIME ZONE 'UTC' as time_start,
to_timestamp(floor((extract('epoch' from res.time_old) / trunc )) * trunc)
AT TIME ZONE 'UTC' + trunc * interval '1 sec' as time_end
from res, t2
GROUP BY res.source, time_start, time_end
--end


--name: UniqueErrors
--meta.cache: true
--input: string date {pattern: datetime}
--connection: System:Datagrok
select et.friendly_name, et.id, count(*)
from event_types et
join events e on e.event_type_id = et.id
join users_sessions s on e.session_id = s.id
join users u on u.id = s.user_id
where @date(e.event_time)
and et.source = 'error'
and et.friendly_name is not null
and et.friendly_name != ''
--and et.is_error = true
group by et.friendly_name, et.id
limit 50;
--end
