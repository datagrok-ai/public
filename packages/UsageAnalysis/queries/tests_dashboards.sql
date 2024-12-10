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
left join test_runs r on r.test_name = t.name and r.build_name = b.name and r.benchmark  and (t.first_run <= r.date_time or t.first_run is null)
where t.name not like '%Unhandled exceptions: Exception' and (not t.name = ': : ')
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
       r.params ->> 'totalWorkers' total_workers,  r.params ->> 'worker' worker
from test_runs r
where r.stress_test and (r.batch_name = COALESCE(NULLIF(@batch_name, ''), (select name  from last_build))) and not r.skipped and not r.test_name like '%Unhandled exceptions: Exception' and not r.test_name = ': : '
group by r.date_time, r.test_name, r.passed, r.duration, r.params ->> 'worker', r.params ->> 'totalWorkers', r.result
order by r.test_name, worker
--end


--name:  LastBuildsBenchmarksCompare
--friendlyName:UA | Tests |  Last Builds Benchmarks Compare
--connection: System:Datagrok
with last_builds as (
    select name  from builds order by build_date desc limit 10
),
benchmarks  as (
  select * from tests 
  where exists
  (select * from test_runs where test_name = name and benchmark)
), 
testruns as (
  select build_name, test_name, duration,
  case when r.passed is null then 'did not run' when r.skipped then 'skipped' when r.passed then 'passed' when not r.passed then 'failed' else 'unknown' end as status
  from test_runs r
  where  r.benchmark = true 
and not r.test_name like '%Unhandled exceptions: Exception'
and r.build_name in (select name from last_builds)
)

select CAST(min(r.duration) AS int)as minDuration, CAST(max(r.duration) AS int)as maxDuration, CAST(avg(r.duration) AS int)as avgDuration,
b.name as test, l.name as build_name, COALESCE(r.status, 'didn''t run') as status
from benchmarks b
join last_builds l on 1=1 
left join testruns as r 
on (r.test_name = b.name and l.name = r.build_name)
group by  b.name, l.name,  r.status
order by l.name desc,  b.name

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
