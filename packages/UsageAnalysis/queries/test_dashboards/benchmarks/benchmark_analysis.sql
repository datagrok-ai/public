--name: Benchmark Analysis
--connection: System:Datagrok
select CAST(min(r.duration) AS int)as duration, r.test_name as test, r.build_name
from test_runs r
where r.benchmark = true
and r.passed = true
and r.skipped = false
and not r.test_name like '%Unhandled exceptions: Exception'
group by r.test_name, r.build_name
order by r.build_name desc, r.test_name