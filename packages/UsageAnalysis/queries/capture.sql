--name: CaptureRules
--friendlyName: Capture Rules
--input: string date {pattern: datetime}
--connection: System:Datagrok
with rules as (
  select r.*,
    case when r.status = 'active' and r.expires_at <= (now() at time zone 'utc') then 'expired' else r.status end as state,
    extract(epoch from coalesce(r.ended_at, r.expires_at) - r.created_at) as span
  from capture_rules r
  where (r.status = 'active' and r.expires_at > (now() at time zone 'utc')) or @date(r.created_at)
)
select 'cap-' || r.number as rule, u.login as author,
  case when r.subject_type = 'everyone' then 'everyone'
    else r.subject_type || ' ' || coalesce(r.subject_name, r.subject_value, '') end
    || case when r.anonymous then ' (anonymous)' else '' end as subject,
  case when r.scope_type is not null then r.scope_type || ' ' || coalesce(r.scope_value, '')
    when r.capture->>'serverLevel' is not null then 'server ' || (r.capture->>'serverLevel')
    else 'all activity' end as scope,
  r.reason,
  case when r.state = 'active' then 'active'
    else r.state || ' after ' || case when r.span >= 86400 then round(r.span / 86400) || ' d'
      when r.span >= 3600 then round(r.span / 3600) || ' h'
      else greatest(0, round(r.span / 60)) || ' min' end end as active,
  r.events_captured as events,
  r.state as status,
  concat_ws(',',
    case when (r.capture->>'clicks')::boolean then 'clicks' end,
    case when (r.capture->>'inputs')::boolean then 'inputs' end,
    case when (r.capture->>'requests')::boolean then 'requests' end,
    case when (r.capture->>'calls')::boolean then 'calls' end,
    case when (r.capture->>'errors')::boolean then 'errors' end,
    case when r.capture->>'serverLevel' is not null then 'server:' || (r.capture->>'serverLevel')
      || coalesce('=' || nullif(array_to_string(array(select jsonb_array_elements_text(
        case when jsonb_typeof(r.capture->'debugFlags') = 'array' then r.capture->'debugFlags' end)), ','), ''), '')
      end) as capture,
  r.anonymous, r.name, r.max_events, r.window_minutes, r.max_sessions,
  r.created_at, r.expires_at, r.ended_at, r.id::text as id, eu.login as stopped_by,
  case when r.state = 'stopped' then (select v.value from events e
    join event_types et on et.id = e.event_type_id
    join event_parameter_values v on v.event_id = e.id
    join event_parameters p on p.id = v.parameter_id and p.name = 'reason'
    where et.source = 'audit' and et.friendly_name = 'capture-rule-stopped'
      and e.event_time between r.ended_at - interval '1 minute' and r.ended_at + interval '1 minute'
      and exists (select 1 from event_parameter_values rv join event_parameters rp on rp.id = rv.parameter_id
        where rv.event_id = e.id and rp.name = 'ruleId' and rv.value = r.id::text)
    limit 1) end as stop_reason
from rules r
left join users u on u.id = r.created_by
left join users eu on eu.id = r.ended_by
order by r.state = 'active' desc, r.number desc
--end
