--name: count all
--tags: log,counts
--connection: System:Datagrok
--output: int eventsCount
--output: int notificationsCount
--output: int usersCount
--output: int projectsCount
--output: int sessionsCount
--output: int groupsCount
--output: int chatsCount
--output: int commentsCount
--output: int connectionsCount
--output: int queriesCount
--output: int jobsCount
--output: int queryRunsCount
--output: int jobRunsCount
--output: int entitiesCount
SELECT
  (SELECT count(*) from events) as eventsCount,
  (SELECT count(*) from notifications) as notificationsCount,
  (SELECT count(*) from users) as usersCount,
  (select count(*) from projects) as projectsCount,
  (SELECT count(*) from users_sessions) as sessionsCount,
  (SELECT count(*) from groups) as groupsCount,
  (SELECT count(*) from chats) as chatsCount,
  (SELECT count(*) from comments) as commentsCount,
  (SELECT count(*) from connections) as connectionsCount,
  (SELECT count(*) from queries) as queriesCount,
  (SELECT count(*) from jobs) as jobsCount,
  (SELECT count(*) from entities) as entitiesCount
