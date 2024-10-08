
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
left join test_runs r on r.test_name = t.name and r.build_name = b.name and
r.date_time = (select max(_r.date_time) from test_runs _r where _r.test_name = r.test_name and _r.build_name = r.build_name and not _r.stress_test)
where not r.stress_test
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
  CAST(min(r.duration) AS int)as duration
from tests t full join last_builds b on 1=1
left join test_runs r on r.test_name = t.name and r.build_name = b.name
where r.benchmark
group by b.name, t.name, r.date_time, r.passed, r.skipped, r.result, t.type
order by b.name, t.name
--end

--name: Stress tests
--friendlyName:UA | Tests | Stress tests
--connection: System:Datagrok
--input: string batch_name {nullable: true; choices:  Query("select distinct batch_name from test_runs where stress_test and NOT batch_name = '' order by 1  desc")}
WITH last_build AS (
    select batch_name as name  from test_runs  where stress_test  order by date_time desc limit 1
)

select r.date_time, 
r.test_name, 
r.passed, 
r.duration, 
r.params ->> 'totalWorkers' total_workers,  r.params ->> 'worker' worker, avg(rs.duration)
from test_runs r
left join test_runs rs on rs.test_name = r.test_name  and not rs.stress_test and rs.passed
where r.stress_test and (r.batch_name = COALESCE(NULLIF(@batch_name, ''), (select name  from last_build)))
group by r.date_time, r.test_name, r.passed, r.duration, r.params ->> 'worker', r.params ->> 'totalWorkers'
order by r.test_name, worker
--end


--name:  LastBuildsBenchmarksCompare
--friendlyName:UA | Tests |  Last Builds Benchmarks Compare
--connection: System:Datagrok
with last_builds as (
    select name  from builds order by build_date desc limit 10
) 

select CAST(min(r.duration) AS int)as duration, r.test_name as test, r.build_name
from test_runs r
where r.benchmark = true
and r.passed = true
and r.skipped = false
and not r.test_name like '%Unhandled exceptions: Exception'
and r.build_name in (select name from last_builds)
group by r.test_name, r.build_name
order by r.build_name desc, r.test_name
--end

--name:  LastBuildsCompare
--friendlyName:UA | Tests | Last Builds Compare
--connection: System:Datagrok
with last_builds as (
    select name  from builds order by build_date desc limit 2
)
select distinct on (b.name, t.name) 
  b.name as build, 
  t.name as test, 
  t.type, 
  t.package, 
  r.date_time,
  case when r.passed is null then 'did not run' when r.skipped then 'skipped' when r.passed then 'passed' when not r.passed then 'failed' else 'unknown' end as status,
  r.result, 
  r.duration, 
  r.benchmark
from tests t 
full join last_builds b  on 1=1
left join test_runs r on r.test_name = t.name and r.build_name = b.name  
order by b.name desc, t.name, r.date_time desc
--end


--name:  LastVersionsCompare
--friendlyName:UA | Tests | Last Versions Compare
--connection: System:Datagrok
with versions_builds as (
  select distinct on ((regexp_matches(name, '\d{4}-\d{2}-\d{2}-(\d+\.\d+\.\d+)'))[1]) 
  name, 
  (regexp_matches(name, '\d{4}-\d{2}-\d{2}-(\d+\.\d+\.\d+)'))[1] as version, 
  build_date 
  from builds 
  order by 
  (regexp_matches(name, '\d{4}-\d{2}-\d{2}-(\d+\.\d+\.\d+)'))[1] desc, 
  build_date desc
),
last_builds as (
  select * 
  from versions_builds 
  order by build_date desc 
  limit 2
) 

select distinct on (b.name, t.name) 
  b.name as build, 
  b.version as version, 
  t.name as test, 
  t.type, 
  t.package, 
  r.date_time,
  case when r.passed is null then 'did not run' when r.skipped then 'skipped' when r.passed then 'passed' when not r.passed then 'failed' else 'unknown' end as status,
  r.result, 
  r.duration, 
  r.benchmark
from tests t full join last_builds b on 1=1
left join test_runs r on r.test_name = t.name and r.build_name = b.name
order by b.name desc, t.name, r.date_time desc
--end
