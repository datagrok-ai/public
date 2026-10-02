--name: TopDisabledErrors
--input: string date { pattern: datetime }
--input: list<string> users
--connection: System:Datagrok
select et.friendly_name, et.id, count(1) from event_types et
join events e on e.event_type_id = et.id
join users_sessions s on e.session_id = s.id
join users u on u.id = s.user_id
where @date(e.event_time)
and (u.login = any(@users) or @users = ARRAY['all'])
and et.source = 'error'
and et.friendly_name is not null
and et.friendly_name != ''
and et.is_error = false
group by et.friendly_name, et.id
limit 50;
--end


--name: TopPackageErrors
--input: string date { pattern: datetime }
--input: list<string> users
--connection: System:Datagrok
select e.error_message, count(1) from event_types et
join entities en on et.id = en.id
join published_packages pp on en.package_id = pp.id
join events e on e.event_type_id = et.id
join users_sessions s on e.session_id = s.id
join users u on u.id = s.user_id
where
e.error_message is not null
and et.friendly_name is not null
and et.friendly_name != ''
and @date(e.event_time)
and (u.login = any(@users) or @users = ARRAY['all'])
group by e.error_message
limit 50;
--end


--name: TopErrorSources
--input: string date { pattern: datetime }
--connection: System:Datagrok
select et.error_source, count(1) as count from event_types et
join events e on e.event_type_id = et.id
where @date(e.event_time)
and et.source = 'error'
group by et.error_source
order by count desc limit 50;
--end


--name: TopErrors
--input: string date {pattern: datetime}
--connection: System:Datagrok
SELECT COALESCE(e.description, t.error_source || ': ' || e.friendly_name || E'\\n' || COALESCE(e.error_stack_trace, ''))
           AS error, COUNT(t.error_stack_trace_hash)
FROM events e
         JOIN event_types t ON e.event_type_id = t.id
WHERE t.source = 'error' AND @date(e.event_time)
GROUP BY error
ORDER BY count DESC LIMIT 15;
--end

--name: ErrorSessions
--friendlyName: Error Sessions
--input: list<string> signatures
--input: list<string> users
--input: string from
--input: string to
--connection: System:Datagrok
select e.session_id::text as session, u.login as "user", min(e.event_time) as first, max(e.event_time) as last,
  count(*)::int as count
from events e
join event_types t on t.id = e.event_type_id
join users_sessions s on s.id = e.session_id
join users u on u.id = s.user_id
where t.source = 'error' and t.error_stack_trace_hash::text = any(@signatures::text[])
  and u.login = any(@users)
  and e.event_time >= @from::timestamp and e.event_time < @to::timestamp
group by 1, 2
order by last desc
limit 50
--end

--name: ErrorReports
--friendlyName: Error Reports
--input: list<string> signatures
--connection: System:Datagrok
select r.number, r.created_on, r.is_auto, r.is_resolved, u.login as reporter,
  left(coalesce(nullif(r.description, ''), r.error_message), 200) as description
from reports r
left join users u on u.id = r.reporter_id
where r.error_stack_trace_hash::text = any(@signatures::text[])
order by r.number desc
limit 50
--end

--name: ErrorAlerts
--friendlyName: Error Alerts
--input: list<string> signatures
--connection: System:Datagrok
select a.kind, a.key, a.status, p.status as problem, a.severity, a.summary, a.opened_at, a.cleared_at, a.resolved_at
from alerts a
left join problems p on p.id = a.problem_id
where (a.kind = 'error-incident' and a.key = any(array(select left(s, 12) from unnest(@signatures::text[]) s)))
  or (a.kind = 'report' and a.key in (select r.id::text from reports r where r.error_stack_trace_hash::text = any(@signatures::text[])))
order by a.opened_at desc
limit 50
--end

--name: ErrorSample
--friendlyName: Error Sample
--input: string signature
--connection: System:Datagrok
select e.event_time as time, coalesce(nullif(e.error_message, ''), e.description) as error,
  coalesce(e.error_stack_trace, t.error_stack_trace) as stack
from events e
join event_types t on t.id = e.event_type_id
where t.source = 'error' and t.error_stack_trace_hash = @signature::uuid
order by e.event_time desc
limit 1
--end
