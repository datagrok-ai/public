--name: UserReports
--input: string date {pattern: datetime}
--connection: System:Datagrok
with cte as (select e.id report_id, e.options ->> 'sequence_id' report_number, e.event_time report_time, e.description,
    e.options ->> 'error_message' as error, e.options ->> 'error_stack_trace' error_stack_trace,
    e.options ->> 'error_stack_trace_hash' error_stack_trace_hash, u.friendly_name reporter
from events e
join event_types t on e.event_type_id = t.id
join users_sessions s on e.session_id = s.id
join users u on u.id = s.user_id
where t.source = 'usage' and t.friendly_name = 'user report posted'
and @date(e.event_time))

select (count(e.id)) same_errors_count, c.report_id, c.report_number, c.report_time, c.description, c.error,
       c.error_stack_trace, c.error_stack_trace_hash, c.reporter
from events e
join event_types t on e.event_type_id = t.id
right join cte c on ((t.error_stack_trace_hash = c.error_stack_trace_hash and c.error_stack_trace_hash is not null) or t.friendly_name = c.error)
group by c.report_id, c.report_number, c.report_time, c.description, c.error,
         c.error_stack_trace, c.error_stack_trace_hash, c.reporter
order by report_time desc
--end
