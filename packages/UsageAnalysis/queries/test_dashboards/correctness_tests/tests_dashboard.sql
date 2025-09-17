--name: TestsDashboard
--friendlyName: UA | Tests | Tests
--connection: System:Datagrok
--input: string instanceFilter = 'dev' {choices: ['', 'dev', 'release', 'public']}
--input: int lastBuildsNum = 5
--input: string versionFilter {nullable: true}
--input: string packageFilter {nullable: true}
--input: string versionFilter {nullable: true}
--input: bool showNotRun = false {optional: true}
--input: bool showBenchmarks = false {optional: true}
--input: bool showNotCiCd = false {optional: true}
WITH last_builds AS (
    select name as build_name, build_date, commit
    from builds b
  	
    -- todo: filter builds with cicd not-stresstest runs
    -- todo: filter only success builds
    where exists (SELECT 1 from test_runs r join tests t on r.test_name = t.name where t.type = 'package'
                                              and r.build_name = b.name
                                              and not r.stress_test
                                              and (@instanceFilter is null or r.instance like 'https://' || @instanceFilter || '.datagrok.ai')
                                              and r.benchmark = @showBenchmarks limit 1
          )
    order by b.build_date desc limit @lastBuildsNum
), last_builds_indexed AS (
  select build_name, ROW_NUMBER() OVER (ORDER BY build_date desc) AS build_index, build_date, commit
  from last_builds b
)

select 
  b.build_name as build,
  b.build_index,
  b.build_date,
  t.name as test,
  t.type, 
  r.date_time,
  r.instance as instance,
  b.commit as build_commit,
  case when r.passed is null then 'did not run' when r.skipped then 'skipped' when r.passed then 'passed' when not r.passed then 'failed' else 'unknown' end as status,
  COALESCE(r.params->>'flaking', 'false')::bool as flaking, 
  r.result,
  r.duration,
  coalesce(nullif(t.owner, ''), p.package_author, '') as owner,
  p.updated_on as last_package_update
from tests t
full join last_builds_indexed b on 1=1
left join test_runs r on 
r.test_name = t.name and 
r.build_name = b.build_name and 
(t.first_run <= r.date_time or t.first_run is null) and
r.date_time = (select max(_r.date_time) from test_runs _r where _r.test_name = r.test_name and _r.build_name = r.build_name and not _r.stress_test) and not r.stress_test and (@showNotCiCd or r.ci_cd)
left join published_packages p on p.name = t.package and p.is_current
where 
t.name not like '%Unhandled exceptions: Exception'
and (not t.name = ': : ')
and not t.type = 'manual'
and (@showNotRun or not r.passed is null)
and r.benchmark = @showBenchmarks
and (@packageFilter is null or @packageFilter = p.name)
and (@versionFilter is null or @versionFilter = b.build_name)
and (@instanceFilter is null or r.instance like 'https://' || @instanceFilter || '.datagrok.ai')
order by b.build_index, t.name