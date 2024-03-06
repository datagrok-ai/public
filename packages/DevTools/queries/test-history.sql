--name: TestHistory
--connection: System:Datagrok
--input: string packageName {nullable :true}
--input: string category {nullable :true}
--input: string test {nullable :true}
select
e.event_time as date,
case when v4.value::bool then 'skipped' when v1.value::bool then 'passed' else 'failed' end as status,
v3.value::int as ms,
v2.value as result,
v5.value as logs
from events e
inner join event_types t on t.id = e.event_type_id and t.source = 'usage' and t.friendly_name like 'test-%'
left join event_parameter_values v1 inner join event_parameters p1 on p1.id = v1.parameter_id and p1.name = 'success' on v1.event_id = e.id
left join event_parameter_values v2 inner join event_parameters p2 on p2.id = v2.parameter_id and p2.name = 'result' on v2.event_id = e.id
left join event_parameter_values v3 inner join event_parameters p3 on p3.id = v3.parameter_id and p3.name = 'ms' on v3.event_id = e.id
left join event_parameter_values v4 inner join event_parameters p4 on p4.id = v4.parameter_id and p4.name = 'skipped' on v4.event_id = e.id
left join event_parameter_values v5 inner join event_parameters p5 on p5.id = v5.parameter_id and p5.name = 'logs' on v5.event_id = e.id
where
@packageName is null or exists(select 1 from event_parameter_values v inner join event_parameters p on p.id = v.parameter_id and p.name = 'packageName' where v.event_id = e.id and v.value = @packageName)
and
@category is null or exists(select 1 from event_parameter_values v inner join event_parameters p on p.id = v.parameter_id and p.name = 'category' where v.event_id = e.id and v.value = @category)
and
@test is null or exists(select 1 from event_parameter_values v inner join event_parameters p on p.id = v.parameter_id and p.name = 'test' where v.event_id = e.id and v.value = @test)
order by e.event_time desc
--end

--name: CategoryHistory
--connection: System:Datagrok
--input: string packageName
--input: string category
select
e.event_time as date, v1.value::bool as success, v2.value as passed, v3.value as skipped, v4.value as failed
from events e
inner join event_types t on t.id = e.event_type_id and t.source = 'usage' and t.friendly_name like 'category-%'
left join event_parameter_values v1 inner join event_parameters p1 on p1.id = v1.parameter_id and p1.name = 'success' on v1.event_id = e.id
left join event_parameter_values v2 inner join event_parameters p2 on p2.id = v2.parameter_id and p2.name = 'passed' on v2.event_id = e.id
left join event_parameter_values v3 inner join event_parameters p3 on p3.id = v3.parameter_id and p3.name = 'skipped' on v3.event_id = e.id
left join event_parameter_values v4 inner join event_parameters p4 on p4.id = v4.parameter_id and p4.name = 'failed' on v4.event_id = e.id
where
@packageName is null or exists(select 1 from event_parameter_values v inner join event_parameters p on p.id = v.parameter_id and p.name = 'packageName' where v.event_id = e.id and v.value = @packageName)
and
@category is null or exists(select 1 from event_parameter_values v inner join event_parameters p on p.id = v.parameter_id and p.name = 'category' where v.event_id = e.id and v.value = @category)
--end

--name: PackageHistory
--connection: System:Datagrok
--input: string packageName
select
e.event_time as date, v1.value::bool as success, v2.value as passed, v3.value as skipped, v4.value as failed
from events e
inner join event_types t on t.id = e.event_type_id and t.source = 'usage' and t.friendly_name like 'package-%'
left join event_parameter_values v1 inner join event_parameters p1 on p1.id = v1.parameter_id and p1.name = 'success' on v1.event_id = e.id
left join event_parameter_values v2 inner join event_parameters p2 on p2.id = v2.parameter_id and p2.name = 'passed' on v2.event_id = e.id
left join event_parameter_values v3 inner join event_parameters p3 on p3.id = v3.parameter_id and p3.name = 'skipped' on v3.event_id = e.id
left join event_parameter_values v4 inner join event_parameters p4 on p4.id = v4.parameter_id and p4.name = 'failed' on v4.event_id = e.id
where
@packageName is null or exists(select 1 from event_parameter_values v inner join event_parameters p on p.id = v.parameter_id and p.name = 'packageName' where v.event_id = e.id and v.value = @packageName)
--end
