--name: StressTestsSummary
--friendlyName: UA | Tests | Stress Tests Summary
--connection: System:Datagrok
--input: int lastBuildsNum = 10
WITH last_builds AS (
    select name, build_date
    from builds b
    where exists (select 1 from stress_tests s where s.build_name = b.name)
    order by b.build_date desc limit @lastBuildsNum
), runs AS (
    select b.name as build, b.build_date, s.total_concurrent_runs as threads,
           s.test_name, s.duration, s.passed
    from last_builds b
        inner join stress_tests s on s.build_name = b.name and not s.skipped
), per_test AS (
    -- each test contributes its own pass fraction (1.0 only if it never failed)
    select build, threads, test_name,
           sum(case when passed then 1 else 0 end)::numeric / count(*) as pass_frac
    from runs group by build, threads, test_name
), pass AS (
    -- average of per-test fractions = sum(pass_frac) / number of tests
    select build, threads, avg(pass_frac) as avg_pass_frac
    from per_test group by build, threads
)
select
    r.build,
    max(r.build_date) as started,
    r.threads,
    count(*) as runs,
    count(distinct r.test_name) as tests,
    min(r.duration) as min_ms,
    cast(round(avg(r.duration)) as int) as avg_ms,
    cast(round(percentile_cont(0.5) within group (order by r.duration)) as int) as median_ms,
    cast(round(percentile_cont(0.95) within group (order by r.duration)) as int) as p95_ms,
    max(r.duration) as max_ms,
    cast(round(100.0 * p.avg_pass_frac) as int) as pass_rate
from runs r
    inner join pass p on p.build = r.build and p.threads = r.threads
group by r.build, r.threads, p.avg_pass_frac
order by started desc, r.threads
--end
