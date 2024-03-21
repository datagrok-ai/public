--name: UserReports
--input: string date {pattern: datetime}
--input: string event_id {nullable: true}
--connection: System:Datagrok
with report_info as (
    select er.report_id report_id, e.event_time report_time, e.description description,
           array_agg(er.event_id) errors_ids, u.friendly_name as reporter from events_reports er
           inner join events e on e.id = er.report_id
           inner join users_sessions s on e.session_id = s.id
           inner join users u on u.id = s.user_id
    group by er.report_id, e.event_time, e.description, u.friendly_name
    having @date(e.event_time)
    )

SELECT r.report_id, r.report_time, r.description,
       case
           when cardinality(r.errors_ids) > 1
               then cardinality(r.errors_ids) || ' ' || 'same errors'
           else '1 error' end as same_errors_count,
       r.reporter,
       COALESCE(e.description, t.error_source || ': ' || e.friendly_name || E'\\n' || COALESCE(e.error_stack_trace, ''))
                              AS error
FROM report_info r
         inner join events e on e.id = r.errors_ids[1]
         inner join event_types t on e.event_type_id = t.id
where r.errors_ids @> array[@event_id]::uuid[]	or @event_id is null
--end