--name: Stress Test Analysis
--connection: System:Datagrok
--input: string batch_name {choices:  Query("select distinct batch_name from test_runs where stress_test and NOT batch_name = '' order by 1  desc")}
select r.date_time, r.test_name, r.passed, r.duration, r.params ->> 'totalWorkers' total_workers,  r.params ->> 'worker' worker, avg(rs.duration) from test_runs r
left join test_runs rs on rs.test_name = r.test_name  and not rs.stress_test and rs.passed
where r.stress_test and r.batch_name = @batch_name
group by r.date_time, r.test_name, r.passed, r.duration, r.params ->> 'worker', r.params ->> 'totalWorkers'
order by r.test_name, worker