--name: # events
--connection: System:Datagrok
--output: int eventsCount
SELECT count(*)as event_count from events ;
--end

--name: # notifications
--connection: System:Datagrok
--output: int notificationsCount
SELECT count(*) as notifications_count from notifications ;
--end

--name: # users
--connection: System:Datagrok
--output: int usersCount
SELECT count(*) from users as users_count;
--end

--name: # projects
--connection: System:Datagrok
--output: int projectsCount
SELECT count(*) as projects_count from projects;
--end

--name: # sessions
--connection: System:Datagrok
--output: int sessionsCount
SELECT count(*) as sessions_count from users_sessions;
--end

--name: # groups
--connection: System:Datagrok
--output: int groupsCount
SELECT count(*)as groups_count from groups ;
--end

--name: # chats
--connection: System:Datagrok
--output: int chatsCount
SELECT count(*) as chats_count from chats;
--end

--name: # comments
--connection: System:Datagrok
--output: int commentsCount
SELECT count(*) as comments_count from comments;
--end

--name: # connections
--connection: System:Datagrok
--output: int connectionsCount
SELECT count(*) as connections_count from connections;
--end

--name: # queries
--connection: System:Datagrok
--output: int queriesCount
SELECT count(*) as queries_count from queries;
--end

--name: # data jobs
--connection: System:Datagrok
--output: int jobsCount
SELECT count(*) as data_jobs_count from jobs;
--end
