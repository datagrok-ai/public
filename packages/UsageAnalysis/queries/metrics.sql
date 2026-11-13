--name: MetricsResetPgStatStatements
--connection: System:DatagrokAdmin
SELECT pg_stat_statements_reset();
--end


--name: MetricsCacheMissTables
--connection: System:Datagrok
--meta.cache: all
--meta.cache.invalidateOn: */5 * * * *
SELECT
  schemaname || '.' || relname AS offender_table,
  ROUND((100.0 * (heap_blks_hit + idx_blks_hit)
    / NULLIF(heap_blks_hit + idx_blks_hit + heap_blks_read + idx_blks_read, 0))::numeric, 1)::float
    AS offender_hit_pct,
  (heap_blks_read + idx_blks_read)::bigint AS offender_disk_reads
FROM pg_statio_user_tables
WHERE heap_blks_read + idx_blks_read >= 1000
  AND (100.0 * (heap_blks_hit + idx_blks_hit)
    / NULLIF(heap_blks_hit + idx_blks_hit + heap_blks_read + idx_blks_read, 0)) < 95
ORDER BY (100.0 * (heap_blks_hit + idx_blks_hit)
  / NULLIF(heap_blks_hit + idx_blks_hit + heap_blks_read + idx_blks_read, 0)) ASC
LIMIT 5;
--end


--name: MetricsTableHealthSummary
--input: int limit = 10
--connection: System:Datagrok
--meta.cache: all
--meta.cache.invalidateOn: */5 * * * *
WITH unhealthy AS (
  SELECT
    schemaname || '.' || relname AS table_name,
    (100.0 * n_dead_tup / NULLIF(n_live_tup + n_dead_tup, 0))::int AS dead_pct,
    GREATEST(last_vacuum, last_autovacuum) AS last_vacuum
  FROM pg_stat_user_tables
  WHERE n_live_tup >= 10000
    AND (100.0 * n_dead_tup / NULLIF(n_live_tup + n_dead_tup, 0)) > 40
)
SELECT
  table_name, dead_pct, last_vacuum,
  COUNT(*) OVER () AS unhealthy_count,
  MAX(dead_pct) OVER () AS max_dead_pct
FROM unhealthy
ORDER BY dead_pct DESC
LIMIT @limit;
--end


--name: MetricsConnectionsOffenders
--input: int limit = 10
--input: int idleXactSec = 60
--input: int activeSec = 30
--connection: System:Datagrok
WITH a AS (
  SELECT
    pid,
    state,
    wait_event_type,
    wait_event,
    application_name,
    usename,
    client_addr::text AS client_addr,
    query,
    EXTRACT(EPOCH FROM (now() - state_change))::int AS state_age_sec,
    EXTRACT(EPOCH FROM (now() - query_start))::int AS query_age_sec,
    cardinality(pg_blocking_pids(pid)) AS blocks
  FROM pg_stat_activity
  WHERE datname = current_database()
    AND backend_type = 'client backend'
    AND pid <> pg_backend_pid()
)
SELECT
  pid, state, wait_event_type, wait_event, application_name, usename, client_addr,
  substring(query, 1, 200) AS query,
  CASE WHEN state = 'idle in transaction' THEN state_age_sec ELSE query_age_sec END AS age_sec,
  blocks
FROM a
WHERE (state = 'idle in transaction' AND state_age_sec >= @idleXactSec)
   OR (state = 'active' AND query_age_sec >= @activeSec)
   OR blocks > 0
ORDER BY blocks DESC, age_sec DESC
LIMIT @limit;
--end


--name: MetricsLargestTables
--input: int limit = 10
--connection: System:Datagrok
--meta.cache: all
--meta.cache.invalidateOn: */5 * * * *
SELECT
  schemaname || '.' || relname AS table_name,
  pg_size_pretty(pg_total_relation_size(relid)) AS total,
  pg_size_pretty(pg_indexes_size(relid)) AS "index",
  n_live_tup AS "#rows",
  pg_total_relation_size(relid) AS total_bytes
FROM pg_stat_user_tables
ORDER BY pg_total_relation_size(relid) DESC
LIMIT @limit;
--end


--name: MetricsTableHealth
--input: int limit = 10
--connection: System:Datagrok
--meta.cache: all
--meta.cache.invalidateOn: */5 * * * *
SELECT
  schemaname || '.' || relname AS table_name,
  COALESCE(round((100.0 * n_dead_tup / NULLIF(n_live_tup + n_dead_tup, 0))::numeric, 0), 0)::int AS dead_pct,
  GREATEST(last_vacuum, last_autovacuum) AS last_vacuum
FROM pg_stat_user_tables
WHERE n_live_tup >= 1000
ORDER BY dead_pct DESC NULLS LAST
LIMIT @limit;
--end


--name: MetricsErrorsCount
--input: string date {pattern: datetime}
--connection: System:Datagrok
--meta.cache: all
--meta.cache.invalidateOn: */5 * * * *
WITH _dates AS (
  SELECT min(event_time) AS min_date, max(event_time) AS max_date FROM events WHERE @date(event_time)
),
dates AS (
  SELECT min_date - (max_date - min_date) AS min_prev_date, min_date, max_date FROM _dates
)
SELECT
  count(*) FILTER (WHERE e.event_time >= d.min_date) AS errors_now,
  count(*) FILTER (WHERE e.event_time < d.min_date) AS errors_prev,
  (SELECT min_date      FROM dates) AS window_start,
  (SELECT max_date      FROM dates) AS window_end,
  (SELECT min_prev_date FROM dates) AS prev_window_start
FROM events e
CROSS JOIN dates d
WHERE e.event_time BETWEEN d.min_prev_date AND d.max_date
  AND e.session_id IS NOT NULL
  AND e.event_type_id IN (SELECT id FROM event_types WHERE source = 'error');
--end


--name: MetricsSessionsCount
--input: string date {pattern: datetime}
--connection: System:Datagrok
--meta.cache: all
--meta.cache.invalidateOn: */5 * * * *
WITH _dates AS (
  SELECT min(event_time) AS min_date, max(event_time) AS max_date FROM events WHERE @date(event_time)
),
dates AS (
  SELECT min_date - (max_date - min_date) AS min_prev_date, min_date, max_date FROM _dates
)
SELECT
  count(DISTINCT s.id) FILTER (WHERE s.started >= d.min_date) AS sessions_now,
  count(DISTINCT s.id) FILTER (WHERE s.started < d.min_date) AS sessions_prev,
  (SELECT min_date      FROM dates) AS window_start,
  (SELECT max_date      FROM dates) AS window_end,
  (SELECT min_prev_date FROM dates) AS prev_window_start
FROM users_sessions s
CROSS JOIN dates d
WHERE s.started BETWEEN d.min_prev_date AND d.max_date;
--end


--name: MetricsRouteFailures
--input: string method
--input: string route
--input: string date {pattern: datetime}
--connection: System:Datagrok
SELECT r.time, r.status, r.ms, u.login AS user, r.session_id,
  err.source AS log_level, err.message AS error, err.stack_trace
FROM http_requests r
JOIN http_routes t ON t.id = r.route_id
LEFT JOIN users_sessions s ON s.id = r.session_id
LEFT JOIN users u ON u.id = s.user_id
LEFT JOIN LATERAL (
  SELECT et.source, coalesce(e.error_message, e.description, e.friendly_name) AS message, e.error_stack_trace AS stack_trace
  FROM events e JOIN event_types et ON et.id = e.event_type_id
  WHERE e.session_id = r.session_id AND et.source IN ('error', 'warning')
    AND e.event_time BETWEEN r.time - interval '5 seconds' AND r.time + interval '5 seconds'
  ORDER BY et.source = 'error' DESC, abs(extract(epoch FROM e.event_time - r.time))
  LIMIT 1
) err ON true
WHERE t.method = @method AND t.route = @route AND r.status >= 400 AND @date(r.time)
ORDER BY r.time DESC
LIMIT 1000;
--end


--name: MetricsRequestLog
--input: string session
--input: string time
--connection: System:Datagrok
SELECT e.event_time AS time, et.source AS level, coalesce(e.error_message, e.description, e.friendly_name) AS message,
  e.error_stack_trace AS stack_trace
FROM events e
JOIN event_types et ON et.id = e.event_type_id
WHERE e.session_id = @session::uuid
  AND e.event_time BETWEEN @time::timestamp - interval '5 seconds' AND @time::timestamp + interval '5 seconds'
ORDER BY e.event_time;
--end
