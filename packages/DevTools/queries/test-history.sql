--name: TestHistory
--connection: System:Datagrok
--input: string packageName {nullable :true}
--input: string category {nullable :true}
--input: string test {nullable :true}
select distinct on (t.name) 
  r.date_time as date, 
  case when r.passed is null then 'did not run' when r.skipped then 'skipped' when r.passed then 'passed' when not r.passed then 'failed' else 'unknown' end as status,
  r.result,
  r.duration as ms, 
  r.params::json->>'logs' as logs
from tests t 
left join test_runs r on r.test_name = t.name  
where t.type='package'
and t.name = CONCAT(@packageName, ': ', @category,': ', @test)
order by t.name, r.date_time desc
--end

--name: CategoryHistory
--connection: System:Datagrok
--input: string packageName
--input: string category
with category_runs as (
select distinct on (b.name, t.name) 
  b.name as build,
  r.date_time as date, 
  case when r.passed is null then 'did not run' when r.skipped then 'skipped' when r.passed then 'passed' when not r.passed then 'failed' else 'unknown' end as status,
  r.result,
  r.duration as ms, 
  r.params::json->>'logs' as logs
from tests t 
full join builds b on 1=1
left join test_runs r on r.test_name = t.name and r.build_name = b.name
where t.type='package'
and t.name like CONCAT(@packageName, ': ', @category,': %')
order by b.name desc, t.name, r.date_time desc)

select  Count(cat_r1.date) as passed, Count(cat_r2.date)as failed, Count(cat_r3.date) as skipped,
case when Min(cat_r1.date) < Min(cat_r2.date) and Min(cat_r1.date) < Min(cat_r3.date) then cat_r1.date when Min(cat_r2.date) < Min(cat_r3.date) then cat_r2.date else cat_r3.date end as date, 
case when Count(cat_r2.date) = 0 then true else false end as success
from builds b 
left join category_runs cat_r1 on cat_r1.build = b.name and cat_r1.status = 'passed'
left join category_runs cat_r2 on cat_r2.build = b.name and cat_r2.status = 'failed'
left join category_runs cat_r3 on cat_r3.build = b.name and cat_r3.status = 'skipped'
group by (b.name, cat_r1.date, cat_r2.date, cat_r3.date)
having not (Count(cat_r3.date)) = 0 or
not (Count(cat_r2.date)) = 0 or
not (Count(cat_r1.date)) = 0
--end

--name: PackageHistory
--connection: System:Datagrok
--input: string packageName
with category_runs as (
select distinct on (b.name, t.name) 
  b.name as build,
  r.date_time as date, 
  case when r.passed is null then 'did not run' when r.skipped then 'skipped' when r.passed then 'passed' when not r.passed then 'failed' else 'unknown' end as status,
  r.result,
  r.duration as ms, 
  r.params::json->>'logs' as logs
from tests t 
full join builds b on 1=1
left join test_runs r on r.test_name = t.name and r.build_name = b.name
where t.type='package'
and t.name like CONCAT(@packageName, ': %')
order by b.name desc, t.name, r.date_time desc)

select  Count(cat_r1.date) as passed, Count(cat_r2.date)as failed, Count(cat_r3.date) as skipped,
case when Min(cat_r1.date) < Min(cat_r2.date) and Min(cat_r1.date) < Min(cat_r3.date) then cat_r1.date when Min(cat_r2.date) < Min(cat_r3.date) then cat_r2.date else cat_r3.date end as date, 
case when Count(cat_r2.date) = 0 then true else false end as success
from builds b 
left join category_runs cat_r1 on cat_r1.build = b.name and cat_r1.status = 'passed'
left join category_runs cat_r2 on cat_r2.build = b.name and cat_r2.status = 'failed'
left join category_runs cat_r3 on cat_r3.build = b.name and cat_r3.status = 'skipped'
group by (b.name, cat_r1.date, cat_r2.date, cat_r3.date)
having not (Count(cat_r3.date)) = 0 or
not (Count(cat_r2.date)) = 0 or
not (Count(cat_r1.date)) = 0
--end
