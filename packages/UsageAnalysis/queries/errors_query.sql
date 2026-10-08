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
where t.source = 'error' and t.error_stack_trace_hash::text = any(@signatures ::text[])
  and u.login = any(@users)
  and e.event_time >= @from ::timestamp and e.event_time < @to ::timestamp
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
where r.error_stack_trace_hash::text = any(@signatures ::text[])
order by r.number desc
limit 50
--end

--name: ErrorAlerts
--friendlyName: Error Alerts
--input: list<string> signatures
--connection: System:Datagrok
select p.kind, p.key, p.alert_status as status, p.status as problem, p.severity, p.summary, p.alerted_at as opened_at,
  case when p.state = 'cleared' then p.cleared_at end as cleared_at, p.resolved_at
from problems p
where p.alert_status is not null
  and ((p.kind = 'error-incident' and p.key = any(array(select left(s, 12) from unnest(@signatures ::text[]) s)))
    or (p.kind = 'report' and p.key in (select r.id::text from reports r where r.error_stack_trace_hash::text = any(@signatures ::text[]))))
order by p.alerted_at desc
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
where t.source = 'error' and t.error_stack_trace_hash = @signature ::uuid
order by e.event_time desc
limit 1
--end
