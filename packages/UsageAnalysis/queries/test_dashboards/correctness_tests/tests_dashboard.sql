--name: TestsDashboard
--friendlyName: UA | Tests | Tests
--connection: System:Datagrok
--input: string instanceFilter = 'dev' {choices: ['', 'dev', 'release', 'public', 'release-ec2']}
--input: int lastBuildsNum = 3
--input: string versionFilter {nullable: true}
--input: string packageFilter {nullable: true}
--input: bool showNotRun = false {optional: true}
--input: bool showBenchmarks = false {optional: true}
--input: bool showNotCiCd = false {optional: true}
WITH last_builds AS (
    SELECT name AS build_name, build_date, commit
    FROM builds b
    WHERE EXISTS (
        SELECT 1 FROM test_runs r
        JOIN tests t ON r.test_name = t.name
        WHERE t.type = 'package'
          AND r.build_name = b.name
          AND NOT r.stress_test
          AND r.benchmark = @showBenchmarks
          AND r.instance LIKE 'https://' || @instanceFilter || '.datagrok.ai'
        LIMIT 1
    )
    ORDER BY b.build_date DESC
    LIMIT @lastBuildsNum
),
last_builds_indexed AS (
    SELECT build_name, ROW_NUMBER() OVER (ORDER BY build_date DESC) AS build_index, build_date, commit
    FROM last_builds
),
latest_runs AS (
    SELECT r.*,
        ROW_NUMBER() OVER (PARTITION BY r.test_name, r.build_name ORDER BY r.date_time DESC) AS rn
    FROM test_runs r
    WHERE NOT r.stress_test
      AND r.benchmark = @showBenchmarks
      AND (@showNotCiCd OR r.ci_cd)
      AND r.instance LIKE 'https://' || @instanceFilter || '.datagrok.ai'
      AND r.build_name IN (SELECT build_name FROM last_builds)
)
SELECT
  b.build_name AS build,
  b.build_index,
  b.build_date,
  t.name AS test,
  t.type,
  r.date_time,
  r.instance AS instance,
  b.commit AS build_commit,
  CASE WHEN r.passed IS NULL THEN 'did not run' WHEN r.skipped THEN 'skipped' WHEN r.passed THEN 'passed' WHEN NOT r.passed THEN 'failed' ELSE 'unknown' END AS status,
  COALESCE(r.params->>'flaking', 'false')::bool AS flaking,
  r.result,
  r.duration,
  COALESCE(NULLIF(t.owner, ''), p.package_author, '') AS owner,
  p.updated_on AS last_package_update
FROM tests t
CROSS JOIN last_builds_indexed b
LEFT JOIN latest_runs r ON r.test_name = t.name AND r.build_name = b.build_name AND r.rn = 1
LEFT JOIN published_packages p ON p.name = t.package AND p.is_current
WHERE t.name NOT LIKE '%Unhandled exceptions: Exception'
  AND t.name != ': : '
  AND t.type != 'manual'
  AND (@showNotRun OR r.passed IS NOT NULL)
  AND (@packageFilter IS NULL OR @packageFilter = p.name)
  AND (@versionFilter IS NULL OR @versionFilter = b.build_name)
ORDER BY b.build_index, t.name
