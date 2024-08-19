--name: UserReports
--input: string date {pattern: datetime}
--connection: System:Datagrok
select id, number, created_on as time,
case when description = 'Auto report' then error_message else description end as description,
error_message as error, error_stack_trace, error_stack_trace_hash, is_resolved, reporter_id as reporter,
    assignee_id as assignee, jira_ticket as jira, array_to_string(labels, ',') as labels
from reports
order by (is_resolved is true) asc, time desc;
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
select id report_id, number, created_on as time,
case when description = 'Auto report' then error_message else description end as description,
error_message as error, error_stack_trace, error_stack_trace_hash, is_resolved, reporter_id as reporter,
    assignee_id as assignee, jira_ticket as jira, array_to_string(labels, ',') as label
from reports
order by (is_resolved is true) asc, time desc limit 20
--end

--name: UserReportsSingle
--input: int reportNumber
--connection: System:Datagrok
select id, number, created_on as time,
case when description = 'Auto report' then error_message else description end as description,
error_message as error, error_stack_trace, error_stack_trace_hash, is_resolved, reporter_id as reporter,
    assignee_id as assignee, jira_ticket as jira, array_to_string(labels, ',') as labels
from reports
where number = @reportNumber;
--end

--name: Reports migration
--connection: System:Datagrok
with cte as (
delete from events e
using event_types t, users_sessions s
where e.event_type_id = t.id and t.source = 'usage' and t.friendly_name = 'user report posted'
and s.id = e.session_id
returning e.id report_id, e.options ->> 'sequence_id' report_number, e.event_time as time,
    e.description,
    e.options ->> 'error_message' as error, e.options ->> 'error_stack_trace' error_stack_trace,
    e.options ->> 'error_stack_trace_hash' error_stack_trace_hash, (e.options ->> 'is_acknowledged')::boolean as is_acknowledged,
    e.options ->> 'assignee_id' as assignee,
    s.user_id reporter, e.options ->> 'jira_ticket_number' as jira, e.options ->> 'label' as label)

INSERT INTO reports (
    id, number,
   created_on,
   type,
   friendly_name,
   name,
   jira_ticket,
   is_auto,
   is_resolved,
   description,
   reporter_id,
   assignee_id,
   error_message,
   error_stack_trace,
   error_stack_trace_hash,
   labels)
SELECT report_id, report_number::int, time, case when label like '%feedback%' then 'feedback' else 'report' end as type,
               null, null, jira, case when description = 'Auto report' then true else false end as is_auto,
               is_acknowledged, description, reporter::uuid, assignee::uuid, error, error_stack_trace, error_stack_trace_hash,
               string_to_array(label, ',')::varchar[] from cte order by time;
--end


--name: ReportDataMigration
--connection: System:Datagrok
--input: string report_id
--input: string id
--input: string screenshot {nullable: true}
--input: string details {nullable: true}
--input: string client_settings {nullable: true}
--input: string server_settings {nullable: true}
--input: string errors {nullable: true}
--input: string client_log {nullable: true}
--input: string server_log {nullable: true}
--input: string console {nullable: true}
--input: string queries_log {nullable: true}
--input: string containers_log {nullable: true}
--input: string images_log {nullable: true}
--input: string services {nullable: true}
with subquery as (
    insert into reports_data (id,
                              screenshot,
                              details,
                              client_settings,
                              server_settings,
                              errors,
                              client_log,
                              server_log,
                              console,
                              queries_log,
                              scripting_log,
                              containers_log,
                              images_log,
                              services)
        values (@id::uuid, @screenshot::jsonb, @details::jsonb, @client_settings::jsonb, @server_settings::jsonb, @errors::jsonb,
        @client_log::jsonb, @server_log::jsonb, @console::jsonb, @queries_log::jsonb, @containers_log::jsonb, @images_log::jsonb, @services::jsonb)
        returning id
)
update reports set data_id = subquery.id
from subquery where reports.id = @report_id
--end
