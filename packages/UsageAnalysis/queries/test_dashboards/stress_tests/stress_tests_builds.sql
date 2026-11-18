--name: StressTestsBuilds
--friendlyName: UA | Tests | Stress Tests Builds
--connection: System:Datagrok
--input: int lastBuildsNum = 10
-- The last trunk builds that reported stress runs, then the sweeps of the latest release branch
-- that has any, with how many rows each holds. The count is the cache key of
-- StressTestsBuildSummary: it comes from the build_name index, so it costs a fraction of the
-- summary itself, and it changes whenever a build is still reporting or was re-run.
WITH stress_builds AS (
    select b.name, b.build_date, coalesce(b.branch, '') as branch,
           b.branch ~ '^release/\d+\.\d+\.\d+$' as is_release
    from builds b
    where exists (select 1 from stress_tests s where s.build_name = b.name)
), release_branch AS (
    select branch from stress_builds where is_release
    order by string_to_array(substring(branch from 9), '.')::int[] desc limit 1
), picked AS (
    (select name, build_date, branch from stress_builds where is_release is not true
     order by build_date desc limit @lastBuildsNum)
    union all
    (select name, build_date, branch from stress_builds
     where branch = (select branch from release_branch)
     order by build_date desc limit 5)
)
select b.name as build, b.build_date as started, b.branch, c.runs
from picked b
    cross join lateral (select cast(count(*) as int) as runs from stress_tests s where s.build_name = b.name) c
order by b.branch like 'release/%', b.build_date desc
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
-- tests_per_s is the level's throughput: each thread's own rate (its tests over its wall time)
-- averaged and multiplied by the thread count. Per thread, because threads of one level start
-- staggered and one build can hold several sweeps, so a level has no single wall-clock span.
-- started is when the build's first sweep began: builds order by it, because several images can
-- share one build date and a sweep can run a day after its image was built.
WITH runs AS (
    select s.total_concurrent_runs as threads, s.test_name, s.duration, s.passed,
           s.task_start, s.task_finish, s.concurrent_run
    from stress_tests s
    where s.build_name = @build and not s.skipped and @runs >= 0
), per_test AS (
    select threads, test_name,
           sum(case when passed then 1 else 0 end)::numeric / count(*) as pass_frac
    from runs group by threads, test_name
), pass AS (
    select threads, avg(pass_frac) as avg_pass_frac
    from per_test group by threads
), tasks AS (
    select threads, count(*) / nullif(extract(epoch from max(task_finish) - task_start), 0) as rate
    from runs group by threads, task_start, concurrent_run
), speed AS (
    select threads, threads * avg(rate) as tests_per_s
    from tasks group by threads
)
select
    cast(@build as varchar) as build,
    (select min(task_start) from runs) as started,
    r.threads,
    count(*) as runs,
    count(distinct r.test_name) as tests,
    min(r.duration) as min_ms,
    cast(round(avg(r.duration)) as int) as avg_ms,
    cast(round(percentile_cont(0.5) within group (order by r.duration)) as int) as median_ms,
    cast(round(percentile_cont(0.95) within group (order by r.duration)) as int) as p95_ms,
    max(r.duration) as max_ms,
    cast(round(100.0 * p.avg_pass_frac) as int) as pass_rate,
    cast(round(cast(sp.tests_per_s as numeric), 1) as double precision) as tests_per_s
from runs r
    inner join pass p on p.threads = r.threads
    inner join speed sp on sp.threads = r.threads
group by r.threads, p.avg_pass_frac, sp.tests_per_s
order by r.threads
--end
