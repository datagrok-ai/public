--name: _groupCompleter
--input: string sub
--connection: System:Datagrok
select id, name from (
select id, name, (case when name = @sub then 10 when name like @sub || '%' then 9 else 1 end) as weight from groups
where personal = false and name ilike '%' || @sub || '%'
order by weight desc, name
limit 50
) a
--end


--name: _userCompleter
--input: string sub
--connection: System:Datagrok
select id, name from (
select id, friendly_name as name, (case when
friendly_name = @sub or login = @sub or first_name = @sub or last_name = @sub then 10
when
friendly_name like @sub || '%'  or login like @sub || '%'  or first_name like @sub || '%'  or last_name like @sub || '%'
then 9 else 1 end) as weight from users
where friendly_name ilike '%' || @sub || '%'
   or login ilike '%' || @sub || '%'
   or first_name ilike '%' || @sub || '%'
   or last_name ilike '%' || @sub || '%'
order by weight desc, name
limit 50
) a
--end


--name: _actionCompleter
--input: string sub
--connection: System:Datagrok
select name from (
select distinct name, (case when name = @sub then 10 when name like @sub || '%' then 9 else 1 end) as weight from event_types
where name ilike '%' || @sub || '%'
order by weight desc, name
limit 50
) a
--end


--name: _queriesCompleter
--input: string sub
--connection: System:Datagrok
select id, name from (
select id, name, (case when name = @sub then 10 when name like @sub || '%' then 9 else 1 end) as weight from queries
where name ilike '%' || @sub || '%'
order by weight desc, name
limit 50
) a
--end


--name: _projectCompleter
--input: string sub
--connection: System:Datagrok
select id, name from (
select id, name, (case when name = @sub then 10 when name like @sub || '%' then 9 else 1 end) as weight  from projects
where name ilike '%' || @sub || '%'
order by weight desc, name
limit 50
) a
--end


--name: @action by @user on @date
--input: string action { suggestions: LogAnalysis:actionCompleter }
--input: string user { suggestions: LogAnalysis:userCompleter}
--input: string date { pattern: datetime }
--connection: System:Datagrok
select t.name, t.source, e.event_time from events e
inner join event_types t on e.event_type_id = t.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where u.id = @user and @date(e.event_time) and t.name = @action
--end


--name: @group members
--tags: log
--input: string group { suggestions: LogAnalysis:groupCompleter }
--connection: System:Datagrok
WITH RECURSIVE children(id, name, is_admin, depth, path, path1, cycle) AS (
        SELECT g.id, g.friendly_name, false, 1, ARRAY[g.id], ''::text, false
        FROM groups g
        where g.id = @group
      UNION ALL
         SELECT g.id, g.friendly_name,  c.is_admin or (case when c.depth = 1 then r.is_admin else false end),
          c.depth + 1, path || g.id, path1 || '/' || g.friendly_name, g.id = ANY(path)
        FROM children c inner join groups_relations r on r.parent_id=c.id
        inner join groups g on r.child_id = g.id
        WHERE NOT cycle)
		select u.name, u.login, c.path1 as path from children c inner join users u on u.group_id = c.id
--end


--name: actions by @user on @date
--input: string user { suggestions: LogAnalysis:userCompleter }
--input: string date { pattern: datetime }
--connection: System:Datagrok
select t.name, t.source, e.event_time from events e
inner join event_types t on e.event_type_id = t.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where u.id = @user and @date(e.event_time)
--end
