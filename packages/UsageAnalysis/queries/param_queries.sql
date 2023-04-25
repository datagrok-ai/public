--name: _groupCompleter
--meta.cache: true
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
--meta.cache: true
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
--meta.cache: true
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
--meta.cache: true
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
--meta.cache: true
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
--meta.cache: true
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
--meta.cache: true
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
--meta.cache: true
--input: string user { suggestions: LogAnalysis:userCompleter }
--input: string date { pattern: datetime }
--connection: System:Datagrok
select t.name, t.source, e.event_time from events e
inner join event_types t on e.event_type_id = t.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where u.id = @user and @date(e.event_time)
--end


--name: @user events by hours
--tags: log
--meta.cache: true
--input: string @user { suggestions: logAnalysis:userCompleter }
--connection: System:Datagrok
select hh as hours_ago, description  from
(select e.description, extract(hour from now() - e.event_time) as hh
from users u
inner join events e
inner join event_parameter_values v on v.event_id = e.id
inner join event_parameters p on p.event_type_id = e.event_type_id and
p.type = 'entity_id' and p.name = 'user' and p.id = v.parameter_id
on  u.id::varchar = v.value
where u.id = @user and e.event_time >= (now() - '1 day'::INTERVAL)) data
order by hours_ago
--end