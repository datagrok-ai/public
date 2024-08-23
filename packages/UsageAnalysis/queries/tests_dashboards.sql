
--name: Tests
--friendlyName: UA | Tests | Tests
--connection: System:Datagrok
WITH last_builds AS (
    select name  from builds order by build_date desc limit 5
)
select
  b.name as build,
  t.name as test,
  t.type, 
  r.date_time,
  case when r.passed is null then 'did not run' when r.skipped then 'skipped' when r.passed then 'passed' when not r.passed then 'failed' else 'unknown' end as status,
  r.result,
  r.duration,
  r.benchmark
from tests t full join last_builds b on 1=1
left join test_runs r on r.test_name = t.name and r.build_name = b.name
order by b.name, t.name
--end

--name: Benchmarks
--friendlyName: UA | Tests | Benchmarks
--connection: System:Datagrok
WITH last_builds AS (
    select name  from builds order by build_date desc limit 5
)
select
  b.name as build,
  t.name as test,
  t.type, 
  r.date_time,
  case when r.passed is null then 'did not run' when r.skipped then 'skipped' when r.passed then 'passed' when not r.passed then 'failed' else 'unknown' end as status,
  r.result,
  r.duration
from tests t full join last_builds b on 1=1
left join test_runs r on r.test_name = t.name and r.build_name = b.name
where r.benchmark
order by b.name, t.name
--end

--name: Stress tests
--friendlyName:UA | Tests | Stress tests
--connection: System:Datagrok
--input: string batch_name {choices:  Query("select distinct batch_name from test_runs where stress_test and NOT batch_name = '' order by 1  desc")}
select r.date_time, 
r.test_name, 
r.passed, 
r.duration, 
r.params ->> 'totalWorkers' total_workers,  r.params ->> 'worker' worker, avg(rs.duration) 
from test_runs r
left join test_runs rs on rs.test_name = r.test_name  and not rs.stress_test and rs.passed
where r.stress_test and r.batch_name = @batch_name
group by r.date_time, r.test_name, r.passed, r.duration, r.params ->> 'worker', r.params ->> 'totalWorkers'
order by r.test_name, worker
--end