--name: GetUsersInGroups
--input: list<string> groups
--connection: System:Datagrok
with recursive selected_groups as (
  select id from groups where name = any(@groups)
  union
  select gr.child_id as id from selected_groups sg
  join groups_relations gr on sg.id = gr.parent_id
)
select u.login from selected_groups sg
join groups g on g.id = sg.id
join users u on g.id = u.group_id
where g.personal = true
and u.login is not null;
--end
