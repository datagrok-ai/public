--name: AiUsageTimeline
--input: string date {pattern: datetime}
--input: list<string> groups
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
),
turns as (
  select e.id, e.event_time, u.friendly_name as user, u.group_id as ugid,
    max(v.value) filter (where p.name = 'feature') as feature,
    max(v.value) filter (where p.name = 'inputCount')::numeric as input_tokens,
    max(v.value) filter (where p.name = 'outputCount')::numeric as output_tokens,
    max(v.value) filter (where p.name = 'cacheReadCount')::numeric as cache_read_tokens,
    max(v.value) filter (where p.name = 'cacheCreationCount')::numeric as cache_creation_tokens,
    max(v.value) filter (where p.name = 'costUsd')::numeric as cost_usd
  from events e
  join event_types et on e.event_type_id = et.id
  join users_sessions s on e.session_id = s.id
  join users u on u.id = s.user_id
  join event_parameter_values v on v.event_id = e.id
  join event_parameters p on p.id = v.parameter_id
  where et.source = 'usage' and et.friendly_name = 'ai-turn'
    and @date(e.event_time)
    and u.group_id in (select id from selected_groups)
  group by e.id, e.event_time, u.friendly_name, u.group_id
),
buckets as (
  select *, case when max(event_time) over () - min(event_time) over () < interval '3 days' then 'hour' else 'day' end as bucket
  from turns
)
select date_trunc(bucket, event_time) as time_start, sum(input_tokens)::bigint as input_tokens,
  sum(output_tokens)::bigint as output_tokens, sum(cost_usd) as cost_usd
from buckets
group by 1
order by 1
--end

--name: AiUsageByUser
--input: string date {pattern: datetime}
--input: list<string> groups
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
),
turns as (
  select e.id, e.event_time, u.friendly_name as user, u.group_id as ugid,
    max(v.value) filter (where p.name = 'feature') as feature,
    max(v.value) filter (where p.name = 'inputCount')::numeric as input_tokens,
    max(v.value) filter (where p.name = 'outputCount')::numeric as output_tokens,
    max(v.value) filter (where p.name = 'cacheReadCount')::numeric as cache_read_tokens,
    max(v.value) filter (where p.name = 'cacheCreationCount')::numeric as cache_creation_tokens,
    max(v.value) filter (where p.name = 'costUsd')::numeric as cost_usd
  from events e
  join event_types et on e.event_type_id = et.id
  join users_sessions s on e.session_id = s.id
  join users u on u.id = s.user_id
  join event_parameter_values v on v.event_id = e.id
  join event_parameters p on p.id = v.parameter_id
  where et.source = 'usage' and et.friendly_name = 'ai-turn'
    and @date(e.event_time)
    and u.group_id in (select id from selected_groups)
  group by e.id, e.event_time, u.friendly_name, u.group_id
)
select "user", ugid, count(*) as turns, sum(input_tokens)::bigint as input_tokens, sum(output_tokens)::bigint as output_tokens,
  sum(cache_read_tokens)::bigint as cache_read_tokens, sum(cache_creation_tokens)::bigint as cache_creation_tokens,
  sum(cost_usd) as cost_usd
from turns
group by "user", ugid
order by cost_usd desc nulls last, input_tokens desc
--end

--name: AiUsageByGroup
--input: string date {pattern: datetime}
--input: list<string> groups
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
),
turns as (
  select e.id, e.event_time, u.friendly_name as user, u.group_id as ugid,
    max(v.value) filter (where p.name = 'feature') as feature,
    max(v.value) filter (where p.name = 'inputCount')::numeric as input_tokens,
    max(v.value) filter (where p.name = 'outputCount')::numeric as output_tokens,
    max(v.value) filter (where p.name = 'cacheReadCount')::numeric as cache_read_tokens,
    max(v.value) filter (where p.name = 'cacheCreationCount')::numeric as cache_creation_tokens,
    max(v.value) filter (where p.name = 'costUsd')::numeric as cost_usd
  from events e
  join event_types et on e.event_type_id = et.id
  join users_sessions s on e.session_id = s.id
  join users u on u.id = s.user_id
  join event_parameter_values v on v.event_id = e.id
  join event_parameters p on p.id = v.parameter_id
  where et.source = 'usage' and et.friendly_name = 'ai-turn'
    and @date(e.event_time)
    and u.group_id in (select id from selected_groups)
  group by e.id, e.event_time, u.friendly_name, u.group_id
),
memberships as (
  select distinct ugid as member, ugid as group_id from turns
  union
  select m.member, gr.parent_id from memberships m
  join groups_relations gr on gr.child_id = m.group_id
)
select g.friendly_name as "group", count(distinct t.user) as users, count(*) as turns, sum(input_tokens)::bigint as input_tokens, sum(output_tokens)::bigint as output_tokens,
  sum(cache_read_tokens)::bigint as cache_read_tokens, sum(cache_creation_tokens)::bigint as cache_creation_tokens,
  sum(cost_usd) as cost_usd
from turns t
join memberships m on m.member = t.ugid
join groups g on g.id = m.group_id
where not g.personal and not g.hidden and g.friendly_name <> 'All users'
group by g.friendly_name
order by cost_usd desc nulls last, input_tokens desc
--end

--name: AiUsageByFeature
--input: string date {pattern: datetime}
--input: list<string> groups
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
),
turns as (
  select e.id, e.event_time, u.friendly_name as user, u.group_id as ugid,
    max(v.value) filter (where p.name = 'feature') as feature,
    max(v.value) filter (where p.name = 'inputCount')::numeric as input_tokens,
    max(v.value) filter (where p.name = 'outputCount')::numeric as output_tokens,
    max(v.value) filter (where p.name = 'cacheReadCount')::numeric as cache_read_tokens,
    max(v.value) filter (where p.name = 'cacheCreationCount')::numeric as cache_creation_tokens,
    max(v.value) filter (where p.name = 'costUsd')::numeric as cost_usd
  from events e
  join event_types et on e.event_type_id = et.id
  join users_sessions s on e.session_id = s.id
  join users u on u.id = s.user_id
  join event_parameter_values v on v.event_id = e.id
  join event_parameters p on p.id = v.parameter_id
  where et.source = 'usage' and et.friendly_name = 'ai-turn'
    and @date(e.event_time)
    and u.group_id in (select id from selected_groups)
  group by e.id, e.event_time, u.friendly_name, u.group_id
)
select feature, count(*) as turns, sum(input_tokens)::bigint as input_tokens, sum(output_tokens)::bigint as output_tokens,
  sum(cache_read_tokens)::bigint as cache_read_tokens, sum(cache_creation_tokens)::bigint as cache_creation_tokens,
  sum(cost_usd) as cost_usd
from turns
group by feature
order by cost_usd desc nulls last, input_tokens desc
--end
