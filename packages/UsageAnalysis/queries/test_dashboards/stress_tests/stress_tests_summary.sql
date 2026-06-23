--name: StressTestsSummary
--friendlyName: UA | Tests | Stress Tests Summary
--connection: System:Datagrok
--input: int lastBuildsNum = 10
WITH last_builds AS (
    select name, build_date
    from builds b
    where exists (select 1 from stress_tests s where s.build_name = b.name)
    order by b.build_date desc limit @lastBuildsNum
)
select
    b.name as build,
    b.build_date as started,
    s.total_concurrent_runs as threads,
    count(*) as runs,
    min(s.duration) as min_ms,
    cast(round(avg(s.duration)) as int) as avg_ms,
    cast(round(percentile_cont(0.5) within group (order by s.duration)) as int) as median_ms,
    cast(round(percentile_cont(0.95) within group (order by s.duration)) as int) as p95_ms,
    max(s.duration) as max_ms,
    cast(round(100.0 * sum(case when s.passed then 1 else 0 end) / count(*)) as int) as pass_rate
from last_builds b
    inner join stress_tests s on s.build_name = b.name and not s.skipped
group by b.name, b.build_date, s.total_concurrent_runs
order by b.build_date desc, threads
--end
