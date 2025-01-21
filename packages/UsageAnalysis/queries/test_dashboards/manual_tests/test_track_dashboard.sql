--name: Manual Tests
--friendlyName: UA | Tests | Manual Tests
--connection: System:Datagrok
--input: int lastBuildsNum = 5
WITH last_builds AS (
    select name, build_date
    from builds b
    -- todo: filter builds with cicd not-stresstest runs
    -- todo: filter only success builds
    where EXISTS((SELECT 1 from test_runs r where r.build_name = b.name and not r.stress_test))
    order by b.build_date desc limit @lastBuildsNum
), last_builds_indexed AS (
  select name, ROW_NUMBER() OVER (ORDER BY build_date DESC) AS build_index,
         build_date
  from last_builds b
)
select
  b.name as build,
  b.build_index,
  b.build_date,
  t.name as test,
  t.type,
  r.date_time,
  case when r.passed is null then 'did not run' when r.skipped then 'skipped' when r.passed then 'passed' when not r.passed then 'failed' else 'unknown' end as status,
  r.result
from tests t full join last_builds_indexed b on 1=1
left join test_runs r on r.test_name = t.name and r.build_name = b.name and
r.date_time = (select max(_r.date_time) from test_runs _r where _r.test_name = r.test_name and _r.build_name = r.build_name)
where 1=1
and t.type = 'manual'
and r.passed is not null
and t.name not like '%Unhandled exceptions: Exception'
and (not t.name = ': : ')

order by b.name, t.name
--end