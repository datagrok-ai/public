--name: TestsDashboard
--friendlyName: UA | Tests | Tests
--connection: System:Datagrok
--input: bool showNotRun = false
--input: int lastBuildsNum = 5
WITH last_builds AS (
    select name  from builds b
    -- todo: filter builds with cicd not-stresstest runs
    -- todo: filter only success builds
    where EXISTS((SELECT 1 from test_runs r where r.build_name = b.name and not r.stress_test))
    order by b.build_date desc limit @lastBuildsNum
)
select
  b.name as build,
  t.name as test,
  t.type, 
  r.date_time,
  case when r.passed is null then 'did not run' when r.skipped then 'skipped' when r.passed then 'passed' when not r.passed then 'failed' else 'unknown' end as status,
  r.result,
  r.duration,
  t.owner
from tests t full join last_builds b on 1=1
left join test_runs r on r.test_name = t.name and r.build_name = b.name and (t.first_run <= r.date_time or t.first_run is null) and
r.date_time = (select max(_r.date_time) from test_runs _r where _r.test_name = r.test_name and _r.build_name = r.build_name and not _r.stress_test) and not r.stress_test
where   t.name not like '%Unhandled exceptions: Exception' and (not t.name = ': : ') and t.name not like 'Unknown:%' and (@showNotRun or not r.passed is null) and r.benchmark = false
order by b.name, t.name
--end