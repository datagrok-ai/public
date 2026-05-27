--name: MetricsPgStatStatementsVersion
--connection: System:Datagrok
SELECT extversion FROM pg_extension WHERE extname = 'pg_stat_statements';
--end


--name: MetricsDbStats
--connection: System:Datagrok
--meta.cache: all
--meta.cache.invalidateOn: 0 */5 * * * *
WITH summary AS (
  SELECT
    pg_database_size(current_database()) AS db_size_bytes,
    pg_size_pretty(pg_database_size(current_database())) AS db_size_pretty,
    COALESCE(round(100.0 * sum(blks_hit)::numeric / NULLIF(sum(blks_hit) + sum(blks_read), 0), 2), 0)::float AS cache_hit_pct,
    max(stats_reset) AS stats_reset
  FROM pg_stat_database
  WHERE datname = current_database()
),
offenders AS (
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
  LIMIT 5
)
SELECT
  s.db_size_bytes, s.db_size_pretty, s.cache_hit_pct, s.stats_reset,
  o.offender_table, o.offender_hit_pct, o.offender_disk_reads
FROM summary s
LEFT JOIN offenders o ON true
ORDER BY o.offender_hit_pct ASC NULLS LAST;
--end


--name: MetricsTableHealthSummary
--input: int limit = 10
--connection: System:Datagrok
--meta.cache: all
--meta.cache.invalidateOn: 0 */5 * * * *
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


--name: MetricsConnections
--connection: System:Datagrok
SELECT
  count(*) AS total,
  count(*) FILTER (WHERE state = 'active') AS active,
  count(*) FILTER (WHERE state = 'idle') AS idle,
  count(*) FILTER (WHERE state = 'idle in transaction') AS idle_in_xact,
  count(*) FILTER (WHERE wait_event_type = 'Lock') AS waiting_on_lock,
  COALESCE(
    EXTRACT(EPOCH FROM (now() - min(state_change) FILTER (WHERE state = 'idle in transaction')))::int,
    0
  ) AS oldest_idle_xact_seconds
FROM pg_stat_activity
WHERE datname = current_database()
  AND backend_type = 'client backend';
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


--name: MetricsTopSlowestQueries
--input: int limit = 10
--connection: System:Datagrok
SELECT
    query,
    calls,
    round(total_exec_time::numeric, 0)::float AS total_ms,
    round(mean_exec_time::numeric, 1)::float AS mean_ms,
    CASE WHEN shared_blks_hit + shared_blks_read > 0
      THEN round((100.0 * shared_blks_hit / (shared_blks_hit + shared_blks_read))::numeric, 0)::int
      ELSE NULL END AS hit_pct
FROM pg_stat_statements
WHERE dbid = (SELECT oid FROM pg_database WHERE datname = current_database())
  AND query NOT ILIKE '%pg_stat_statements%'
  AND query NOT ILIKE '%pg_stat_user_tables%'
  AND query NOT ILIKE '%pg_stat_activity%'
  AND query NOT ILIKE '%pg_stat_database%'
  AND query NOT ILIKE '%pg_database_size%'
  AND query NOT ILIKE 'begin%'
  AND query NOT ILIKE 'commit%'
  AND query NOT ILIKE 'rollback%'
  AND query <> 'select $1'
ORDER BY total_exec_time DESC
LIMIT @limit;
--end


--name: MetricsTopMostCalledQueries
--input: int limit = 10
--connection: System:Datagrok
SELECT
    query,
    calls,
    round(total_exec_time::numeric, 0)::float AS total_ms,
    round(mean_exec_time::numeric, 1)::float AS mean_ms,
    CASE WHEN shared_blks_hit + shared_blks_read > 0
      THEN round((100.0 * shared_blks_hit / (shared_blks_hit + shared_blks_read))::numeric, 0)::int
      ELSE NULL END AS hit_pct
FROM pg_stat_statements
WHERE dbid = (SELECT oid FROM pg_database WHERE datname = current_database())
  AND query NOT ILIKE '%pg_stat_statements%'
  AND query NOT ILIKE '%pg_stat_user_tables%'
  AND query NOT ILIKE '%pg_stat_activity%'
  AND query NOT ILIKE '%pg_stat_database%'
  AND query NOT ILIKE '%pg_database_size%'
  AND query NOT ILIKE 'begin%'
  AND query NOT ILIKE 'commit%'
  AND query NOT ILIKE 'rollback%'
  AND query <> 'select $1'
ORDER BY calls DESC
LIMIT @limit;
--end


--name: MetricsWorstCacheHitQueries
--input: int limit = 10
--connection: System:Datagrok
SELECT
    query,
    calls,
    round(total_exec_time::numeric, 0)::float AS total_ms,
    round(mean_exec_time::numeric, 1)::float AS mean_ms,
    round((100.0 * shared_blks_hit / NULLIF(shared_blks_hit + shared_blks_read, 0))::numeric, 0)::int AS hit_pct
FROM pg_stat_statements
WHERE dbid = (SELECT oid FROM pg_database WHERE datname = current_database())
  AND query NOT ILIKE '%pg_stat_statements%'
  AND query NOT ILIKE '%pg_stat_user_tables%'
  AND query NOT ILIKE '%pg_stat_activity%'
  AND query NOT ILIKE '%pg_stat_database%'
  AND query NOT ILIKE '%pg_database_size%'
  AND query NOT ILIKE 'begin%'
  AND query NOT ILIKE 'commit%'
  AND query NOT ILIKE 'rollback%'
  AND query <> 'select $1'
  AND calls >= 1000
  AND shared_blks_hit + shared_blks_read > 0
ORDER BY (1.0 * shared_blks_hit / NULLIF(shared_blks_hit + shared_blks_read, 0)) ASC NULLS LAST
LIMIT @limit;
--end


-- PG 12 / PSS 1.7- variants: total_time / mean_time instead of total_exec_time / mean_exec_time.

--name: MetricsTopSlowestQueriesPg12
--input: int limit = 10
--connection: System:Datagrok
SELECT
    query,
    calls,
    round(total_time::numeric, 0)::float AS total_ms,
    round(mean_time::numeric, 1)::float AS mean_ms,
    CASE WHEN shared_blks_hit + shared_blks_read > 0
      THEN round((100.0 * shared_blks_hit / (shared_blks_hit + shared_blks_read))::numeric, 0)::int
      ELSE NULL END AS hit_pct
FROM pg_stat_statements
WHERE dbid = (SELECT oid FROM pg_database WHERE datname = current_database())
  AND query NOT ILIKE '%pg_stat_statements%'
  AND query NOT ILIKE '%pg_stat_user_tables%'
  AND query NOT ILIKE '%pg_stat_activity%'
  AND query NOT ILIKE '%pg_stat_database%'
  AND query NOT ILIKE '%pg_database_size%'
  AND query NOT ILIKE 'begin%'
  AND query NOT ILIKE 'commit%'
  AND query NOT ILIKE 'rollback%'
  AND query <> 'select $1'
ORDER BY total_time DESC
LIMIT @limit;
--end


--name: MetricsTopMostCalledQueriesPg12
--input: int limit = 10
--connection: System:Datagrok
SELECT
    query,
    calls,
    round(total_time::numeric, 0)::float AS total_ms,
    round(mean_time::numeric, 1)::float AS mean_ms,
    CASE WHEN shared_blks_hit + shared_blks_read > 0
      THEN round((100.0 * shared_blks_hit / (shared_blks_hit + shared_blks_read))::numeric, 0)::int
      ELSE NULL END AS hit_pct
FROM pg_stat_statements
WHERE dbid = (SELECT oid FROM pg_database WHERE datname = current_database())
  AND query NOT ILIKE '%pg_stat_statements%'
  AND query NOT ILIKE '%pg_stat_user_tables%'
  AND query NOT ILIKE '%pg_stat_activity%'
  AND query NOT ILIKE '%pg_stat_database%'
  AND query NOT ILIKE '%pg_database_size%'
  AND query NOT ILIKE 'begin%'
  AND query NOT ILIKE 'commit%'
  AND query NOT ILIKE 'rollback%'
  AND query <> 'select $1'
ORDER BY calls DESC
LIMIT @limit;
--end


--name: MetricsWorstCacheHitQueriesPg12
--input: int limit = 10
--connection: System:Datagrok
SELECT
    query,
    calls,
    round(total_time::numeric, 0)::float AS total_ms,
    round(mean_time::numeric, 1)::float AS mean_ms,
    round((100.0 * shared_blks_hit / NULLIF(shared_blks_hit + shared_blks_read, 0))::numeric, 0)::int AS hit_pct
FROM pg_stat_statements
WHERE dbid = (SELECT oid FROM pg_database WHERE datname = current_database())
  AND query NOT ILIKE '%pg_stat_statements%'
  AND query NOT ILIKE '%pg_stat_user_tables%'
  AND query NOT ILIKE '%pg_stat_activity%'
  AND query NOT ILIKE '%pg_stat_database%'
  AND query NOT ILIKE '%pg_database_size%'
  AND query NOT ILIKE 'begin%'
  AND query NOT ILIKE 'commit%'
  AND query NOT ILIKE 'rollback%'
  AND query <> 'select $1'
  AND calls >= 1000
  AND shared_blks_hit + shared_blks_read > 0
ORDER BY (1.0 * shared_blks_hit / NULLIF(shared_blks_hit + shared_blks_read, 0)) ASC NULLS LAST
LIMIT @limit;
--end


--name: MetricsLargestTables
--input: int limit = 10
--connection: System:Datagrok
--meta.cache: all
--meta.cache.invalidateOn: 0 */5 * * * *
SELECT
  schemaname || '.' || relname AS table_name,
  pg_size_pretty(pg_total_relation_size(relid)) AS total,
  pg_size_pretty(pg_indexes_size(relid)) AS "index",
  pg_total_relation_size(relid) AS total_bytes
FROM pg_stat_user_tables
ORDER BY pg_total_relation_size(relid) DESC
LIMIT @limit;
--end


--name: MetricsTableHealth
--input: int limit = 10
--connection: System:Datagrok
--meta.cache: all
--meta.cache.invalidateOn: 0 */5 * * * *
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
--meta.cache.invalidateOn: 0 */5 * * * *
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
--meta.cache.invalidateOn: 0 */5 * * * *
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


--name: MetricsLatency
--input: string date {pattern: datetime}
--connection: System:Datagrok
WITH _dates AS (
  SELECT min(event_time) AS min_date, max(event_time) AS max_date FROM events WHERE @date(event_time)
),
dates AS (
  SELECT min_date - (max_date - min_date) AS min_prev_date, min_date, max_date FROM _dates
),
http_events AS (
  SELECT e.event_time, epv.value::int AS ms
  FROM events e
  CROSS JOIN dates d
  JOIN event_types et ON et.id = e.event_type_id
    AND et.source = 'usage' AND et.friendly_name = 'http-request'
  JOIN event_parameters ep ON ep.event_type_id = et.id AND ep.name = 'ms'
  JOIN event_parameter_values epv ON epv.event_id = e.id AND epv.parameter_id = ep.id
  WHERE e.event_time BETWEEN d.min_prev_date AND d.max_date
)
SELECT
  COALESCE(round(percentile_cont(0.50) WITHIN GROUP (ORDER BY ms)
    FILTER (WHERE event_time >= (SELECT min_date FROM dates)))::int, 0) AS p50_now,
  COALESCE(round(percentile_cont(0.95) WITHIN GROUP (ORDER BY ms)
    FILTER (WHERE event_time >= (SELECT min_date FROM dates)))::int, 0) AS p95_now,
  COALESCE(round(percentile_cont(0.99) WITHIN GROUP (ORDER BY ms)
    FILTER (WHERE event_time >= (SELECT min_date FROM dates)))::int, 0) AS p99_now,
  COALESCE(round(percentile_cont(0.95) WITHIN GROUP (ORDER BY ms)
    FILTER (WHERE event_time <  (SELECT min_date FROM dates)))::int, 0) AS p95_prev,
  count(*) FILTER (WHERE event_time >= (SELECT min_date FROM dates)) AS count_now,
  count(*) FILTER (WHERE event_time <  (SELECT min_date FROM dates)) AS count_prev,
  (SELECT min_date      FROM dates) AS window_start,
  (SELECT max_date      FROM dates) AS window_end,
  (SELECT min_prev_date FROM dates) AS prev_window_start
FROM http_events;
--end


--name: MetricsHttpRoutes
--input: string date {pattern: datetime}
--input: int limit = 10
--connection: System:Datagrok
WITH http_events AS (
  SELECT
    e.id,
    max(epv.value) FILTER (WHERE ep.name = 'method')          AS method,
    max(epv.value) FILTER (WHERE ep.name = 'route')           AS route,
    max(epv.value::int) FILTER (WHERE ep.name = 'status')     AS status,
    max(epv.value::int) FILTER (WHERE ep.name = 'ms')         AS ms
  FROM events e
  JOIN event_types et ON et.id = e.event_type_id
    AND et.source = 'usage' AND et.friendly_name = 'http-request'
  JOIN event_parameter_values epv ON epv.event_id = e.id
  JOIN event_parameters ep ON ep.id = epv.parameter_id
  WHERE @date(e.event_time)
  GROUP BY e.id
)
SELECT
  method || ' ' || route AS route,
  count(*)::int AS count,
  round(percentile_cont(0.50) WITHIN GROUP (ORDER BY ms))::int AS p50,
  round(percentile_cont(0.95) WITHIN GROUP (ORDER BY ms))::int AS p95,
  round(percentile_cont(0.99) WITHIN GROUP (ORDER BY ms))::int AS p99,
  round(100.0 * count(*) FILTER (WHERE status >= 400) / NULLIF(count(*), 0), 1)::float AS err_pct
FROM http_events
WHERE route IS NOT NULL AND ms IS NOT NULL
GROUP BY method, route
ORDER BY p95 DESC NULLS LAST
LIMIT @limit;
--end
