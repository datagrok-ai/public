--name:UniqueUsersOverview
--input: string date {pattern: datetime}
--input: list groups
--meta.cache: true
--meta.invalidate: 0 * * * *
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
)
select date(e.event_time) as date, count(distinct u.id)
	from events e
	inner join users_sessions s on e.session_id = s.id
	inner join users u on u.id = s.user_id
	inner join selected_groups sg on u.group_id = sg.id
	WHERE @date(e.event_time)
group by date(e.event_time)
--end

--name:UniqueUsersStats
--input: string date {pattern: datetime}
--input: list groups
--meta.cache: true
--meta.invalidate: 0 * * * *
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups
  where id = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
)
select u.name as name, u.id as uid, count(e.id) as count
	from events e
	inner join users_sessions s on e.session_id = s.id
	inner join users u on u.id = s.user_id
	inner join selected_groups sg on u.group_id = sg.id
	WHERE @date(e.event_time)
group by u.name, u.id
--end
