--name: FunctionsUsage
--input: string date {pattern: datetime}
--input: list<string> groups
--input: list<string> packages
--input: list<string> tags
--meta.cache: all
--meta.cache.invalidateOn: 0 0 * * *
--connection: System:Datagrok
--test: FunctionsUsage(date='today', ['1ab8b38d-9c4e-4b1e-81c3-ae2bde3e12c5'], ['all'], ['any'])
with recursive selected_groups as (
    select id from groups
    where id = any(@groups)
union
select gr.child_id as id from selected_groups sg
                                  join groups_relations gr on sg.id = gr.parent_id
    ),
res AS (
select DISTINCT e.id as id_, et.name as function,
coalesce(pp.name, p1n.name, p1f.name, 'Core') as package,
u.friendly_name as user, e.event_time as time_old, u.id as uid,
u.group_id as ugid, coalesce(pp.package_id, p1n.id, p1f.id, '00000000-0000-0000-0000-000000000000') as pid
from events e
inner join event_parameter_values evv on evv.event_id = e.id and evv.value_uuid is not null
inner join funcs et on evv.value_uuid = et.id
inner join entities en on et.id = en.id
left join published_packages pp on en.package_id = pp.id
left join project_relations pr ON pr.entity_id = en.id
left join projects proj ON proj.id = pr.project_id
and proj.is_root = true
and proj.is_package = true
left join packages p1n on p1n.name = proj.name
left join packages p1f on p1f.friendly_name = proj.name and p1n.id is null
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
inner join selected_groups sg on u.group_id = sg.id
where @date(e.event_time)
and (@tags = ARRAY['any'] or exists (select 1 from tags t where t.entity_id = en.id and t.tag = any(@tags)))
),
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
select res.function, res.package, res.user, count(*) AS count,
to_timestamp(floor((extract('epoch' from res.time_old) / trunc )) * trunc)
AT TIME ZONE 'UTC' as time_start,
to_timestamp(floor((extract('epoch' from res.time_old) / trunc )) * trunc)
AT TIME ZONE 'UTC' + trunc * interval '1 sec' as time_end,
res.uid, res.ugid, coalesce(res.pid, '00000000-0000-0000-0000-000000000000') as pid
from res, t2
where (res.package = any(@packages) or @packages = ARRAY['all'])
GROUP BY res.function, res.package, res.user, time_start, time_end,
    res.uid, res.ugid, res.pid
--end


--name: FunctionsContextPane
--input: int time_start
--input: int time_end
--input: list<string> users
--input: list<string> packages
--input: list<string> functions
--meta.cache: all
--meta.cache.invalidateOn: 0 0 * * *
--connection: System:Datagrok
with res AS (
    select DISTINCT e.id as id_, coalesce(pp.name, p1n.name, p1f.name, 'Core') as package,
    e.friendly_name as run, et.name as function, e.event_time as time, e.id as rid,
    coalesce(pp.package_id, p1n.id, p1f.id, '00000000-0000-0000-0000-000000000000') as pid
    from events e
    inner join event_parameter_values evv on evv.event_id = e.id and evv.value_uuid is not null
    inner join funcs et on evv.value_uuid = et.id
    inner join entities en on et.id = en.id
    left join published_packages pp on en.package_id = pp.id
    left join project_relations pr ON pr.entity_id = en.id
    left join projects proj ON proj.id = pr.project_id
    and proj.is_root = true
    and proj.is_package = true
    left join packages p1n on p1n.name = proj.name
    left join packages p1f on p1f.friendly_name = proj.name and p1n.id is null
    inner join users_sessions s on e.session_id = s.id
    inner join users u on u.id = s.user_id
    where e.event_time between to_timestamp(@time_start)
    and to_timestamp(@time_end)
    and u.id = any(@users)
    and et.name = any(@functions)
    )
select res.package, res.run, res.function, res.time, res.rid, res.pid
from res
where res.pid = any(@packages)
--end


--name: FunctionsExecTime
--input: string func
--connection: System:Datagrok
with res as (
    select EXTRACT(EPOCH FROM (fc.ended - fc.started)) AS time,
    fp.name as input, fp.type as type, COALESCE(fcv.value, fcv.value_uuid::text, fcv.value_array::text, fcv.data_frame_value_name, fcv.value_string) as value
    from func_calls fc
    join funcs f on f.id = fc.func_id
    left join func_params fp on fp.func_id = f.id
    right join func_call_values fcv on fcv.func_param_id = fp.id and fcv.func_call_id = fc.id
    where f.name = @func
    )
select res.time, json_object_agg(res.input, res.value)::text as input
from res
GROUP BY res.time
--end

--name: EntitiesTags
--connection: System:Datagrok
select DISTINCT tag from tags
--end