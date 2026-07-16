--name: ReleaseTests
--friendlyName: Release | Tests
--connection: System:Datagrok
--input: string instanceFilter = "dev" {choices: ["dev", "release", "public", "release-ec2"]}
--input: int lastBuildsNum = 5
WITH recent AS (
  SELECT b.name AS build_name, b.build_date, b.commit
  FROM builds b
  WHERE EXISTS (
    SELECT 1 FROM test_runs r
    JOIN tests t ON r.test_name = t.name
    WHERE t.type = 'package'
      AND r.build_name = b.name
      AND NOT r.stress_test
      AND NOT r.benchmark
      AND r.instance LIKE 'https://' || @instanceFilter || '.datagrok.ai'
  )
  ORDER BY b.build_date DESC
  LIMIT @lastBuildsNum
),
indexed_builds AS (
  SELECT build_name, build_date, commit,
         CAST(row_number() OVER (ORDER BY build_date) AS int) AS build_index
  FROM recent
),
latest_runs AS (
  -- One row per (test, build): the newest CI/CD run on the selected instance. Collapsing retries/reruns
  -- here (rn = 1) is what keeps a retried-then-passed test from being counted as failed, and drops
  -- non-CI (local/ad-hoc) runs — matching the correctness TestsDashboard.
  SELECT r.*,
         ROW_NUMBER() OVER (PARTITION BY r.test_name, r.build_name ORDER BY r.date_time DESC) AS rn
  FROM test_runs r
  WHERE NOT r.stress_test AND NOT r.benchmark AND r.ci_cd
    AND r.instance LIKE 'https://' || @instanceFilter || '.datagrok.ai'
    AND r.build_name IN (SELECT build_name FROM recent)
),
active_tests AS (
  -- Tests that ran on this instance in the last 14 days, with their most recent CI run date.
  -- Broadened beyond the build window so tests that recently STOPPED running are still surfaced
  -- (for the "not run for >= N days" staleness alert); last_run drives that check client-side.
  SELECT r.test_name, MAX(r.date_time) AS last_run
  FROM test_runs r
  WHERE NOT r.stress_test AND NOT r.benchmark AND r.ci_cd AND r.passed IS NOT NULL
    AND r.instance LIKE 'https://' || @instanceFilter || '.datagrok.ai'
    AND r.date_time >= now() - interval '14 days'
  GROUP BY r.test_name
),
prev_version AS (
  -- Highest release minor strictly below the current window's version (e.g. 1.28 -> 1.27), found across
  -- ALL instances: a previous release ran on its own instance at the time, not the current one.
  SELECT v.minor FROM (
    SELECT DISTINCT (regexp_match(b.name, '(\d+\.\d+)\.'))[1] AS minor
    FROM builds b WHERE (regexp_match(b.name, '(\d+\.\d+)\.'))[1] IS NOT NULL
  ) v
  WHERE string_to_array(v.minor, '.')::int[] <
        string_to_array((SELECT (regexp_match(build_name, '(\d+\.\d+)\.'))[1]
                         FROM indexed_builds ORDER BY build_index DESC LIMIT 1), '.')::int[]
  ORDER BY string_to_array(v.minor, '.')::int[] DESC
  LIMIT 1
),
prev_builds AS (
  SELECT b.name FROM builds b
  WHERE (regexp_match(b.name, '(\d+\.\d+)\.'))[1] = (SELECT minor FROM prev_version)
    AND EXISTS (SELECT 1 FROM test_runs r WHERE r.build_name = b.name
                AND r.passed AND NOT r.stress_test AND NOT r.benchmark)
  ORDER BY b.build_date DESC LIMIT 5
),
prev_avg AS (
  -- Average passed-run duration per test across the previous release's last builds (the baseline).
  SELECT r.test_name, CAST(AVG(r.duration) AS int) AS prev_ms
  FROM test_runs r JOIN prev_builds pb ON pb.name = r.build_name
  WHERE r.passed AND NOT r.stress_test AND NOT r.benchmark
  GROUP BY r.test_name
)
SELECT
  COALESCE(p.name, t.package) AS package,
  t.name AS test,
  t.owner,
  b.build_index,
  b.build_name AS build,
  b.commit AS build_commit,
  b.build_date,
  r.instance,
  act.last_run,
  pa.prev_ms AS prev_release_ms,
  CASE WHEN r.passed IS NULL THEN 'did not run'
       WHEN r.skipped THEN 'skipped'
       WHEN r.passed THEN 'passed'
       WHEN NOT r.passed THEN 'failed'
       ELSE 'unknown' END AS status,
  r.result,
  CAST(r.duration AS int) AS duration,
  COALESCE(r.params->>'flaking', 'false')::bool AS flaking
FROM tests t
CROSS JOIN indexed_builds b
LEFT JOIN latest_runs r
  ON r.test_name = t.name AND r.build_name = b.build_name AND r.rn = 1
LEFT JOIN active_tests act ON act.test_name = t.name
LEFT JOIN prev_avg pa ON pa.test_name = t.name
LEFT JOIN published_packages p ON p.name = t.package AND p.is_current
WHERE t.type <> 'manual'
  AND t.name IN (SELECT test_name FROM active_tests)
  AND t.name NOT LIKE '%Unhandled exceptions: Exception'
  AND t.name <> ': : '
ORDER BY b.build_index, package, t.name
--end

--name: ReleaseManualTests
--friendlyName: Release | Manual Tests
--connection: System:Datagrok
--input: int lastBatchesNum = 5
--meta.cache: all
--meta.cache.invalidateOn: 0 0 * * *
WITH recent_batches AS (
  SELECT r.batch_name, MAX(r.date_time) AS latest_batch_run
  FROM test_runs r
  JOIN tests t ON t.name = r.test_name
  WHERE t.type = 'manual' AND r.batch_name IS NOT NULL AND r.batch_name <> '' AND r.passed IS NOT NULL
  GROUP BY r.batch_name
  ORDER BY latest_batch_run DESC
  LIMIT @lastBatchesNum
),
indexed AS (
  SELECT batch_name, latest_batch_run,
         CAST(row_number() OVER (ORDER BY latest_batch_run DESC) AS int) AS batch_index
  FROM recent_batches
),
active_manual AS (
  SELECT DISTINCT test_name FROM test_runs
  WHERE batch_name IN (SELECT batch_name FROM indexed) AND passed IS NOT NULL
)
SELECT
  b.batch_name,
  b.batch_index,
  t.name AS test,
  CASE WHEN r.passed IS NULL THEN 'did not run'
       WHEN r.skipped THEN 'skipped'
       WHEN r.passed THEN 'passed'
       WHEN NOT r.passed THEN 'failed'
       ELSE 'unknown' END AS status,
  r.result
FROM tests t
CROSS JOIN indexed b
LEFT JOIN test_runs r ON r.test_name = t.name AND r.batch_name = b.batch_name
  AND r.date_time = (SELECT MAX(_r.date_time) FROM test_runs _r
                     WHERE _r.test_name = t.name AND _r.batch_name = b.batch_name)
WHERE t.type = 'manual'
  AND t.name IN (SELECT test_name FROM active_manual)
  AND t.name NOT LIKE '%Unhandled exceptions: Exception'
  AND t.name <> ': : '
ORDER BY b.batch_index, t.name
--end
