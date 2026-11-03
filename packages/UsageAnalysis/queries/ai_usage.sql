--name: AiUsage
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
turn_type as (
  select id from event_types where source = 'usage' and friendly_name = 'ai-turn'
),
params as (
  select p.id, p.name from event_parameters p
  where p.event_type_id in (select id from turn_type)
    and p.name in ('feature', 'model', 'status', 'chatId', 'turnId', 'inputCount', 'outputCount', 'cacheReadCount',
      'cacheCreationCount', 'costUsd', 'durationMs')
),
turn_events as (
  select e.id, e.event_time, s.user_id, u.friendly_name as user_name
  from events e
  join users_sessions s on s.id = e.session_id
  join users u on u.id = s.user_id
  where e.event_type_id in (select id from turn_type)
    and @date(e.event_time)
    and u.group_id in (select id from selected_groups)
),
turns as (
  select e.id, e.event_time, e.user_id, e.user_name, x.*
  from turn_events e
  cross join lateral (
    select
      max(v.value) filter (where p.name = 'feature') as feature,
      max(v.value) filter (where p.name = 'model') as model,
      max(v.value) filter (where p.name = 'status') as status,
      max(v.value) filter (where p.name = 'chatId') as chat_id,
      max(v.value) filter (where p.name = 'turnId') as turn_id,
      max(v.value) filter (where p.name = 'inputCount')::int as input_tokens,
      max(v.value) filter (where p.name = 'outputCount')::int as output_tokens,
      max(v.value) filter (where p.name = 'cacheReadCount')::int as cache_read_tokens,
      max(v.value) filter (where p.name = 'cacheCreationCount')::int as cache_write_tokens,
      max(v.value) filter (where p.name = 'costUsd')::numeric as cost_usd,
      max(v.value) filter (where p.name = 'durationMs')::int as duration_ms
    from event_parameter_values v
    join params p on p.id = v.parameter_id
    where v.event_id = e.id
  ) x
),
memberships as (
  select u.id as user_id, u.group_id from users u
  where u.id in (select user_id from turns)
  union
  select m.user_id, gr.parent_id from memberships m
  join groups_relations gr on gr.child_id = m.group_id
),
user_groups as (
  select m.user_id, string_agg(g.friendly_name, ', ' order by g.friendly_name) as groups
  from memberships m
  join groups g on g.id = m.group_id
  where not g.personal and not g.hidden and g.friendly_name <> 'All users'
  group by m.user_id
)
select t.event_time as time, date_trunc('day', t.event_time) as day, t.user_name as "user", ug.groups,
  t.feature, coalesce(t.model, 'unknown') as model, coalesce(t.status, 'ok') as status, t.chat_id, t.turn_id,
  t.input_tokens, t.output_tokens, t.cache_read_tokens, t.cache_write_tokens,
  coalesce(t.input_tokens, 0) + coalesce(t.output_tokens, 0) + coalesce(t.cache_read_tokens, 0) + coalesce(t.cache_write_tokens, 0) as total_tokens,
  t.cost_usd, t.duration_ms
from turns t
left join user_groups ug on ug.user_id = t.user_id
order by t.event_time desc
--end

--name: AiTurnContent
--input: string turnId
--connection: System:Datagrok
with content_params as (
  select p.id, p.name from event_parameters p
  join event_types et on et.id = p.event_type_id
  where et.source = 'usage' and et.friendly_name = 'ai-turn-content'
)
select
  max(v.value) filter (where p.name = 'prompt') as prompt,
  max(v.value) filter (where p.name = 'response') as response,
  max(v.value) filter (where p.name = 'toolCalls') as tool_calls
from event_parameter_values v
join content_params p on p.id = v.parameter_id
where p.name in ('prompt', 'response', 'toolCalls')
  and v.event_id in (
    select t.event_id from event_parameter_values t
    where t.parameter_id in (select id from content_params where name = 'turnId') and t.value::text = @turnId::text)
--end
