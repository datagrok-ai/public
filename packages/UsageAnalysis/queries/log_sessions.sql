--name: Log Sessions
--connection: System:Datagrok

select u.login, u.picture, us.started, us.ip from users u, users_sessions us
where us.user_id = u.id