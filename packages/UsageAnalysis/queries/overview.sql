--name:UniqueUsersOverview
--input: string date {pattern: datetime}
--input: list<string> groups
--input: list<string> packages
--meta.cache: all
--meta.cache.invalidateOn: 0 0 * * *
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
),
trunc_val as (
  select case
    when mx - mn >= INTERVAL '6 month' then 604800
    when mx - mn >= INTERVAL '70 day' then 172800
    when mx - mn >= INTERVAL '10 day' then 43200
    when mx - mn >= INTERVAL '2 day' then 14400
    when mx - mn >= INTERVAL '3 hour' then 1800
    else 600
  end as trunc
  from (
    select min(event_time) as mn, max(event_time) as mx
    from events where @date(event_time)
  ) t
),
package_event_ids as (
  select epv.event_id
  from published_packages pp
  join event_parameter_values epv on epv.value_uuid = pp.id
  join event_parameters ep on epv.parameter_id = ep.id and ep.name = 'package'
  where pp.name = any(@packages)
  union
  select epv.event_id
  from published_packages pp
  join packages p on p.id = pp.package_id
  join entities e1 on e1.package_id = pp.id
  join event_parameter_values epv on epv.value_uuid = e1.id and epv.value <> 'null'
  join event_parameters ep on epv.parameter_id = ep.id and ep.type = 'entity_id'
  where pp.name = any(@packages)
)
select date, count(*) from (
  select distinct
    to_timestamp(floor(extract('epoch' from e.event_time) / tv.trunc) * tv.trunc)
      AT TIME ZONE 'UTC' as date,
    u.group_id as ugid
  from events e
  join users_sessions s on e.session_id = s.id
  join users u on u.id = s.user_id
  join selected_groups sg on u.group_id = sg.id
  cross join trunc_val tv
  where @date(e.event_time)
    and (@packages = ARRAY['all'] or e.id in (select event_id from package_event_ids))
) sub
group by date
--end

--name: PackagesUsageOverview
--input: string date {pattern: datetime}
--input: list<string> groups
--input: list<string> packages
--meta.cache: all
--meta.cache.invalidateOn: 0 0 * * *
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
),
package_events as (
  select distinct on (sub.event_id) sub.event_id, sub.package_name, sub.package_id
  from (
    select epv.event_id, pp.name as package_name, pp.package_id, 1 as priority
    from published_packages pp
    join event_parameter_values epv on epv.value_uuid = pp.id
    join event_parameters ep on epv.parameter_id = ep.id and ep.name = 'package'
    union all
    select epv.event_id, pp.name, pp.package_id, 2 as priority
    from published_packages pp
    join packages p on p.id = pp.package_id
    join entities e1 on e1.package_id = pp.id
    join event_parameter_values epv on epv.value_uuid = e1.id and epv.value <> 'null'
    join event_parameters ep on epv.parameter_id = ep.id and ep.type = 'entity_id'
  ) sub
  order by sub.event_id, sub.priority
)
select
  coalesce(pe.package_name, 'Core') as package,
  u.friendly_name as user, u.friendly_name as name,
  count(*) as count, u.id as uid, u.group_id as ugid,
  coalesce(pe.package_id, '00000000-0000-0000-0000-000000000000') as pid,
  max(e.event_time) as time_end, min(e.event_time) as time_start
from events e
join users_sessions s on e.session_id = s.id
join users u on u.id = s.user_id
join selected_groups sg on u.group_id = sg.id
left join package_events pe on pe.event_id = e.id
where @date(e.event_time)
  and (coalesce(pe.package_name, 'Core') = any(@packages) or @packages = ARRAY['all'])
group by coalesce(pe.package_name, 'Core'),
  coalesce(pe.package_id, '00000000-0000-0000-0000-000000000000'),
  u.friendly_name, u.id, u.group_id
--end
