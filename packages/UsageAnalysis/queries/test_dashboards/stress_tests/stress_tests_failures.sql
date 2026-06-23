--name: StressTestsFailures
--friendlyName: UA | Tests | Stress Tests Failures
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
    count(*) as runs,
    sum(case when not s.passed then 1 else 0 end) as fails,
    cast(round(100.0 * sum(case when s.passed then 1 else 0 end) / count(*)) as int) as pass_rate,
    max(case when not s.passed then s.result end) as sample_error
from stress_tests s
where s.build_name = (select name from selected) and not s.skipped
    and s.test_name != ': : '
group by s.test_name, s.total_concurrent_runs
having sum(case when not s.passed then 1 else 0 end) > 0
order by s.test_name, s.total_concurrent_runs
--end
