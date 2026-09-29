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
    else r.subject_type || ' ' || coalesce(r.subject_name, r.subject_value, '') end as subject,
  case when r.scope_type is not null then r.scope_type || ' ' || coalesce(r.scope_value, '')
    when r.capture->>'serverLevel' is not null then 'server ' || (r.capture->>'serverLevel')
    else 'all activity' end as scope,
  r.reason,
  case when r.span >= 86400 then round(r.span / 86400) || ' d'
    when r.span >= 3600 then round(r.span / 3600) || ' h'
    else greatest(0, round(r.span / 60)) || ' min' end
    || case when r.state <> 'active' then ' (' || r.state || ')' else '' end as active,
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
  r.created_at, r.expires_at, r.ended_at, r.id::text as id
from rules r
left join users u on u.id = r.created_by
order by r.number desc
--end
