--name: # events
--connection: System:DatagrokAdmin
--output: int eventsCount
SELECT count(*)as event_count from events ;
--end

--name: # notifications
--connection: System:DatagrokAdmin
--output: int notificationsCount
SELECT count(*) as notifications_count from notifications ;
--end

--name: # users
--connection: System:DatagrokAdmin
--output: int usersCount
SELECT count(*) from users as users_count;
--end

--name: # projects
--connection: System:DatagrokAdmin
--output: int projectsCount
SELECT count(*) as projects_count from projects;
--end

--name: # sessions
--connection: System:DatagrokAdmin
--output: int sessionsCount
SELECT count(*) as sessions_count from users_sessions;
--end

--name: # groups
--connection: System:DatagrokAdmin
--output: int groupsCount
SELECT count(*)as groups_count from groups ;
--end

--name: # chats
--connection: System:DatagrokAdmin
--output: int chatsCount
SELECT count(*) as chats_count from chats;
--end

--name: # comments
--connection: System:DatagrokAdmin
--output: int commentsCount
SELECT count(*) as comments_count from comments;
--end

--name: # connections
--connection: System:DatagrokAdmin
--output: int connectionsCount
SELECT count(*) as connections_count from connections;
--end

--name: # queries
--connection: System:DatagrokAdmin
--output: int queriesCount
SELECT count(*) as queries_count from queries;
--end

--name: # data jobs
--connection: System:DatagrokAdmin
--output: int jobsCount
SELECT count(*) as data_jobs_count from jobs;
--end
