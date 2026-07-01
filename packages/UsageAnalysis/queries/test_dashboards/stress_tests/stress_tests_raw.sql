--name: StressTestsRaw
--friendlyName: UA | Tests | Stress Tests Raw
--connection: System:Datagrok
--input: string build {nullable: true; choices: Query("select b.name from builds b where exists (select 1 from stress_tests s where s.build_name = b.name) order by b.build_date desc")}
WITH selected AS (
    select coalesce(nullif(@build, ''),
        (select b.name from builds b where exists (select 1 from stress_tests s where s.build_name = b.name)
         order by b.build_date desc limit 1)) as name
)
select
    s.test_name as test,
    s.total_concurrent_runs as threads,
    s.concurrent_run,
    s.duration as ms,
    s.passed,
    s.result as error,
    s.date_time as started
from stress_tests s
where s.build_name = (select name from selected) and not s.skipped
    and s.test_name not like '%Unhandled exceptions: Exception' and s.test_name != ': : '
order by s.test_name, s.total_concurrent_runs, s.concurrent_run
--end
