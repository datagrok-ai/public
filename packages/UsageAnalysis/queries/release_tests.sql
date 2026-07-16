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
  SELECT DISTINCT test_name FROM latest_runs WHERE rn = 1 AND passed IS NOT NULL
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
