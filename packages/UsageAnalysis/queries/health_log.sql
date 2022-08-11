--name: Save Health Log
--input: dataframe health
--connection: System:TestTrack
CREATE TABLE IF NOT EXISTS health_log (
  server varchar(64),
  key varchar(64),
  status varchar(64),
  enabled bool,
  time timestamp without time zone,
  description varchar(1024)
  );

insert into health_log(server, key, description, status, enabled, time)
select '$GROK_SERVER_API_HOST', key, description, status, enabled, time from health h
where not exists(select * from health_log l where l.key = h.key and l.time = h.time and l.server = '$GROK_SERVER_API_HOST')
--end

--name: Health Log Summary
--input: string eventTime = "today" {pattern: datetime}
--connection: System:TestTrack
select key, (select count(*) from health_log _l where _l.key = l.key and status ='Running') running_events,
(select count(*) from health_log _l where _l.key = l.key and status ='Failed') failed_events
from health_log  l
where @eventTime(l.time) and server = '$GROK_SERVER_API_HOST'
group by key
--end