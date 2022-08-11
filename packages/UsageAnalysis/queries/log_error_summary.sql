--name: Log Error Summary
--input: string eventTime = "today" {pattern: datetime}
--connection: System:Datagrok

select v.value, count(*) from events e
inner join event_parameter_values v on v.event_id = e.id
inner join event_types t on t.id = e.event_type_id
inner join event_parameters p on p.event_type_id = t.id and p.id = v.parameter_id and p.name = 'data'

 where t.source ='audit' and t.name ='error' and @eventTime(e.event_time)
 group by v.value
--end