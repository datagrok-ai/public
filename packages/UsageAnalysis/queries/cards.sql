--name: UniqueUsersCount
--input: string date {pattern: datetime}
--input: list groups
--input: list packages
--output: int count1
--output: int count2
--meta.invalidate: 0 * * * *
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
),
_dates as (
select min(event_time) as min_date, max(event_time) as max_date from events WHERE @date(event_time)
),
dates as (
select min_date - (max_date - min_date) as min_prev_date, min_date, max_date from _dates
),
res AS (
select u.id as uid, case when e.event_time < min_date then 1 else 2 end as period
from events e
inner join event_types et on e.event_type_id = et.id
left join entities en on et.id = en.id
left join published_packages pp on en.package_id = pp.id
left join event_parameter_values epv inner join event_parameters ep on epv.parameter_id = ep.id and ep.name = 'package'
on epv.event_id = e.id
left join published_packages pp1 on pp1.id::text = epv.value
left join event_parameter_values epv1 inner join event_parameters ep1 on epv1.parameter_id = ep1.id and ep1.type = 'entity_id'
inner join entities e1 on epv1.value != 'null' and e1.id = epv1.value::uuid
inner join published_packages pp2 inner join packages p2 on p2.id = pp2.package_id on e1.package_id = pp2.id
on epv1.event_id = e.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
inner join dates d on e.event_time between d.min_prev_date and d.max_date
inner join selected_groups sg on u.group_id = sg.id
where (coalesce(pp.name, pp1.name, pp2.name, 'Core') = any(@packages) or @packages = ARRAY['all'])
)
select (select count(distinct res.uid) as count1 from res where period = 1),
(select count(distinct res.uid) as count2 from res where period = 2)
--end


--name: NewUsersCount
--input: string date {pattern: datetime}
--input: list groups
--output: int count1
--output: int count2
--meta.invalidate: 0 * * * *
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
),
_dates as (
select min(u.joined) as min_date, max(u.joined) as max_date from users u WHERE @date(u.joined)
),
dates as (
select min_date - (max_date - min_date) as min_prev_date, min_date, max_date from _dates
),
res AS (
select u.id as uid, case when u.joined < min_date then 1 else 2 end as period
from users u
inner join dates d on u.joined between d.min_prev_date and d.max_date
inner join selected_groups sg on u.group_id = sg.id
)
select (select count(distinct res.uid) as count1 from res where period = 1),
(select count(distinct res.uid) as count2 from res where period = 2)
--end


--name: SessionsCount
--input: string date {pattern: datetime}
--input: list groups
--input: list packages
--output: int count1
--output: int count2
--meta.invalidate: 0 * * * *
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
),
_dates as (
select min(event_time) as min_date, max(event_time) as max_date from events WHERE @date(event_time)
),
dates as (
select min_date - (max_date - min_date) as min_prev_date, min_date, max_date from _dates
),
res AS (
select s.id as sid, case when e.event_time < min_date then 1 else 2 end as period
from events e
inner join event_types et on e.event_type_id = et.id
left join entities en on et.id = en.id
left join published_packages pp on en.package_id = pp.id
left join event_parameter_values epv inner join event_parameters ep on epv.parameter_id = ep.id and ep.name = 'package'
on epv.event_id = e.id
left join published_packages pp1 on pp1.id::text = epv.value
left join event_parameter_values epv1 inner join event_parameters ep1 on epv1.parameter_id = ep1.id and ep1.type = 'entity_id'
inner join entities e1 on epv1.value != 'null' and e1.id = epv1.value::uuid
inner join published_packages pp2 inner join packages p2 on p2.id = pp2.package_id on e1.package_id = pp2.id
on epv1.event_id = e.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
inner join dates d on e.event_time between d.min_prev_date and d.max_date
inner join selected_groups sg on u.group_id = sg.id
where (coalesce(pp.name, pp1.name, pp2.name, 'Core') = any(@packages) or @packages = ARRAY['all'])
)
select (select count(distinct res.sid) as count1 from res where period = 1),
(select count(distinct res.sid) as count2 from res where period = 2)
--end


--name: ViewsCount
--input: string date {pattern: datetime}
--input: list groups
--output: int count1
--output: int count2
--meta.invalidate: 0 * * * *
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
),
_dates as (
select min(q.created_on) as min_date, max(q.created_on) as max_date
from view_layouts q WHERE @date(q.created_on)
),
dates as (
select min_date - (max_date - min_date) as min_prev_date, min_date, max_date from _dates
),
res AS (
select q.id as qid, case when q.created_on < min_date then 1 else 2 end as period
from view_layouts q
inner join users u on u.id = q.author_id
inner join dates d on q.created_on between d.min_prev_date and d.max_date
inner join selected_groups sg on u.group_id = sg.id
)
select (select count(distinct res.qid) as count1 from res where period = 1),
(select count(distinct res.qid) as count2 from res where period = 2)
--end


--name: ConnectionsCount
--input: string date {pattern: datetime}
--input: list groups
--output: int count1
--output: int count2
--meta.invalidate: 0 * * * *
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
),
_dates as (
select min(c.created_on) as min_date, max(c.created_on) as max_date
from connections c WHERE @date(c.created_on)
),
dates as (
select min_date - (max_date - min_date) as min_prev_date, min_date, max_date from _dates
),
res AS (
select c.id as cid, case when c.created_on < min_date then 1 else 2 end as period
from connections c
inner join users u on u.id = c.author_id
inner join dates d on c.created_on between d.min_prev_date and d.max_date
inner join selected_groups sg on u.group_id = sg.id
)
select (select count(distinct res.cid) as count1 from res where period = 1),
(select count(distinct res.cid) as count2 from res where period = 2)
--end


--name: QueriesCount
--input: string date {pattern: datetime}
--input: list groups
--output: int count1
--output: int count2
--meta.invalidate: 0 * * * *
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
),
_dates as (
select min(q.created_on) as min_date, max(q.created_on) as max_date
from queries q WHERE @date(q.created_on)
),
dates as (
select min_date - (max_date - min_date) as min_prev_date, min_date, max_date from _dates
),
res AS (
select q.id as qid, case when q.created_on < min_date then 1 else 2 end as period
from queries q
inner join users u on u.id = q.author_id
inner join dates d on q.created_on between d.min_prev_date and d.max_date
inner join selected_groups sg on u.group_id = sg.id
)
select (select count(distinct res.qid) as count1 from res where period = 1),
(select count(distinct res.qid) as count2 from res where period = 2)
--end
