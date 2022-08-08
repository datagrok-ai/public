--name: latest_scenario_results
--connection: System:TestTrack
--input: string date = today {pattern: datetime}
select ts.*, ta.date, ta1.result, ta1.env
from test_scenario ts
left join
(
  select scenario_id, max(date) as date
  from test_activity
  where @date(date)
  group by scenario_id
) as ta on ta.scenario_id = ts.id
left join test_activity ta1 on ta.date = ta1.date
--end


--name: scenario_activity
--connection: System:TestTrack
--input: string id
select date, result, error_desc as error
from test_activity
where scenario_id = @id
order by date desc
--end


--name: insert_activity
--connection: System:TestTrack
--input: string scenarioId
--input: string result
--input: string desc = ""
--input: string env {choices: ["dev", "test", "prod"]}
insert into test_activity (scenario_id, date, result, error_desc, env)
values (@scenarioId, now(), @result, @desc, @env)
--end


--name: post_scenario
--connection: System:TestTrack
--input: string scenarioId
--input: string name
--input: string path
--input: string type {choices: ["auto", "manual", "grok"]}
insert into test_scenario(id, file_name, file_path, type)
select @scenarioId, @name, @path, @type where not exists(select 1 from test_scenario where id = @scenarioId)
--end


--name: find_scenario
--connection: System:TestTrack
--input: string name
--input: string path
select * from test_scenario where file_name = @name and file_path = @path
--end


--name: manual activity by @date
--connection: System:TestTrack
--input: string date { pattern: datetime }
select ta.date, ts.file_name, ta.result, ta.error_desc as error
from test_activity ta
join test_scenario ts on ts.id = ta.scenario_id
where @date(ta.date) and ts.type = 'manual'
order by ta.date desc
--end
