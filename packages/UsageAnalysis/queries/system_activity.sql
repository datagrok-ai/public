--name: SystemActivity
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
res as (
  select e.id, e.event_time, et.friendly_name as event, e.description, s.user_id as session_user_id
  from events e
  join event_types et on e.event_type_id = et.id
  left join users_sessions s on e.session_id = s.id
  where et.source = 'audit'
    and et.friendly_name in ('server-started', 'user-logged-in', 'user-logged-out', 'user-login-failed',
      'user-impersonated', 'impersonation-failed', 'admin-session-started', 'admin-session-ended',
      'dev-key-generated', 'settings-changed', 'log-settings-changed')
    and @date(e.event_time)
),
params as (
  select v.event_id,
    (array_agg(v.value_uuid) filter (where p.name = 'user'))[1] as user_id,
    (array_agg(v.value_uuid) filter (where p.name = 'subject'))[1] as subject_id,
    (array_agg(v.value) filter (where p.name = 'login'))[1] as login,
    (array_agg(v.value) filter (where p.name = 'type'))[1] as type,
    (array_agg(v.value) filter (where p.name = 'reason'))[1] as reason,
    (array_agg(v.value) filter (where p.name = 'changed'))[1] as changed,
    (array_agg(v.value) filter (where p.name = 'version'))[1] as version
  from event_parameter_values v
  join event_parameters p on p.id = v.parameter_id
  where v.event_id in (select id from res)
  group by v.event_id
)
select res.event_time, res.event,
  coalesce(u.friendly_name, p.login) as user,
  replace(replace(res.description, '@user', coalesce(u.friendly_name, 'Unknown user')),
    '@subject', coalesce(t.friendly_name, p.subject_id::text, '')) as description,
  coalesce(p.reason, p.changed, p.version, p.type) as details,
  u.group_id as ugid, res.id
from res
left join params p on p.event_id = res.id
left join users u on u.id = coalesce(p.user_id, res.session_user_id)
left join users t on t.id = p.subject_id
where u.id is null or u.group_id in (select id from selected_groups)
order by res.event_time desc
--end


--name: SystemActivitySummary
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
res as (
  select et.friendly_name as event, e.event_time as time_old
  from events e
  join event_types et on e.event_type_id = et.id
  left join users_sessions s on e.session_id = s.id
  left join users u on u.id = s.user_id
  where et.source = 'audit'
    and et.friendly_name in ('server-started', 'user-logged-in', 'user-logged-out', 'user-login-failed',
      'user-impersonated', 'impersonation-failed', 'admin-session-started', 'admin-session-ended',
      'dev-key-generated', 'settings-changed', 'log-settings-changed')
    and @date(e.event_time)
    and (u.id is null or u.group_id in (select id from selected_groups))
),
t1 as (
  select (max(res.time_old) - min(res.time_old)) as inter from res
),
t2 as (
  select case when inter >= interval '6 month' then 864000
    when inter >= interval '70 day' then 216000
    when inter >= interval '10 day' then 86400
    when inter >= interval '2 day' then 14400
    when inter >= interval '3 hour' then 3600
    else 600 end as trunc
  from t1
)
select res.event, count(*),
  to_timestamp(floor((extract('epoch' from res.time_old) / trunc)) * trunc) at time zone 'UTC' as time_start
from res, t2
group by res.event, time_start
order by time_start
--end
