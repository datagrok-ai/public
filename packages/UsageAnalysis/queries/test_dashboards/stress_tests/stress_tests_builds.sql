--name: StressTestsBuilds
--friendlyName: UA | Tests | Stress Tests Builds
--connection: System:Datagrok
--input: int lastBuildsNum = 10
-- The last builds that reported stress runs, with how many rows each holds. The count is the cache
-- key of StressTestsBuildSummary: it comes from the build_name index, so it costs a fraction of the
-- summary itself, and it changes whenever a build is still reporting or was re-run.
select b.name as build, b.build_date as started, c.runs
from (
    select name, build_date from builds b
    where exists (select 1 from stress_tests s where s.build_name = b.name)
    order by b.build_date desc limit @lastBuildsNum
) b
    cross join lateral (select cast(count(*) as int) as runs from stress_tests s where s.build_name = b.name) c
order by b.build_date desc
--end

--name: StressTestsBuildSummary
--friendlyName: UA | Tests | Stress Tests Build Summary
--connection: System:Datagrok
--input: string build
--input: int runs
--meta.cache: all
-- StressTestsSummary for one build. A build holds ~160k rows, so summarizing the last 20 from
-- scratch takes minutes; cached per build, only a new or changed build is ever computed. @runs
-- (from StressTestsBuilds) is there for the cache key alone: same build, same row count, same answer.
WITH runs AS (
    select s.total_concurrent_runs as threads, s.test_name, s.duration, s.passed
    from stress_tests s
    where s.build_name = @build and not s.skipped and @runs >= 0
), per_test AS (
    select threads, test_name,
           sum(case when passed then 1 else 0 end)::numeric / count(*) as pass_frac
    from runs group by threads, test_name
), pass AS (
    select threads, avg(pass_frac) as avg_pass_frac
    from per_test group by threads
)
select
    cast(@build as varchar) as build,
    (select b.build_date from builds b where b.name = @build) as started,
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
    inner join pass p on p.threads = r.threads
group by r.threads, p.avg_pass_frac
order by r.threads
--end
