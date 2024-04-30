--name: UserReports
--input: string date {pattern: datetime}
--connection: System:Datagrok
with cte as (select e.id report_id, e.options ->> 'sequence_id' report_number, e.event_time as time, e.description,
    e.options ->> 'error_message' as error, e.options ->> 'error_stack_trace' error_stack_trace,
    e.options ->> 'error_stack_trace_hash' error_stack_trace_hash, (e.options ->> 'is_acknowledged')::boolean as is_acknowledged,
    u2.friendly_name as assignee,
    u.friendly_name reporter, e.options ->> 'jira_ticket_number' as jira, e.options ->> 'label' as label
from events e
    join event_types t on e.event_type_id = t.id
    join users_sessions s on e.session_id = s.id
    join users u on u.id = s.user_id
    left join users u2 on u2.id = (e.options ->> 'assignee_id')::uuid
where t.source = 'usage' and t.friendly_name = 'user report posted'
group by e.id, u.friendly_name, u2.friendly_name)

select (count(e.id)) errors, c.report_id, c.report_number, c.time, c.description, c.error,
       c.error_stack_trace, c.error_stack_trace_hash, c.is_acknowledged, c.reporter, c.assignee, c.jira, c.label
from cte c
         left join event_types t on ((t.error_stack_trace_hash = c.error_stack_trace_hash and c.error_stack_trace_hash is not null) or t.friendly_name = c.error)
         left join events e on e.event_type_id = t.id
group by c.report_id, c.report_number, c.time, c.description, c.error,
         c.error_stack_trace, c.error_stack_trace_hash, c.is_acknowledged, c.reporter, c.assignee, c.jira, c.label
order by (is_acknowledged is true) asc, errors desc, time desc
--end

--name: ReportSameErrors
--input: string stackTraceHash
--input: string errorMessage
--connection: System:Datagrok
select e.event_time, e.source, e.description, e.error_stack_trace, u.id as user_id, u.friendly_name as user_name, e.id as event_id
FROM events e
JOIN users_sessions s on e.session_id = s.id
JOIN users u on u.id = s.user_id
JOIN event_types t ON e.event_type_id = t.id
WHERE t.error_stack_trace_hash = @stackTraceHash or t.friendly_name = @errorMessage;
--end


--name: ReportsTop20
--input: string packageOwnerId
--connection: System:Datagrok
with cte as (select e.id report_id, e.options ->> 'sequence_id' report_number, e.event_time as time, e.description,
    e.options ->> 'error_message' as error, e.options ->> 'error_stack_trace' error_stack_trace,
    e.options ->> 'error_stack_trace_hash' error_stack_trace_hash, (e.options ->> 'is_acknowledged')::boolean as is_acknowledged,
    u2.friendly_name as assignee,
    u.friendly_name reporter, e.options ->> 'jira_ticket_number' as jira, e.options ->> 'label' as label, e.options ->> 'package_owner_id' package_owner
from events e
    join event_types t on e.event_type_id = t.id
    join users_sessions s on e.session_id = s.id
    join users u on u.id = s.user_id
    left join users u2 on u2.id = (e.options ->> 'assignee_id')::uuid
where t.source = 'usage' and t.friendly_name = 'user report posted'
group by e.id, u.friendly_name, u2.friendly_name)

select (count(e.id)) errors, c.report_id, c.report_number, c.time, c.description, c.is_acknowledged,
       c.reporter, c.label, c.assignee, c.jira, c.package_owner, c.error
from cte c
         left join event_types t on ((t.error_stack_trace_hash = c.error_stack_trace_hash and c.error_stack_trace_hash is not null) or t.friendly_name = c.error)
         left join events e on e.event_type_id = t.id
group by c.report_id, c.report_number, c.time, c.description, c.error,
         c.error_stack_trace, c.error_stack_trace_hash, c.is_acknowledged, c.reporter, c.assignee, c.jira, c.label, c.package_owner
order by (is_acknowledged is true) asc, (case when c.package_owner = @packageOwnerId then 1 else 0 end) desc, errors desc, time desc
    limit 20
--end