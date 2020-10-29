--name: total users count
--connection: datagrok
select count(*) from users;
--end


--name: new users
--input: string interv { pattern: datetime }
--connection: datagrok
select count(*) as count from users
where @interv(joined);
--end


--name: new users last month week day
--connection: datagrok
with
month_count as (select count(*) as month_count from users where joined > current_timestamp - interval '30 days'),
week_count as (select count(*) as week_count from users where joined > current_timestamp - interval '7 days'),
day_count as (select count(*) as day_count from users where joined > current_timestamp - interval '1 day')
select * from month_count, week_count, day_count;
--end


--name: new events and errors last month week day
--connection: datagrok
with
month_events as (select count(*) as month_events from events where event_time > current_timestamp - interval '30 days'),
month_errors as (select count(*) as month_errors from events where error_stack_trace is not null and event_time > current_timestamp - interval '30 days'),

week_events as (select count(*) as week_events from events where event_time > current_timestamp - interval '7 days'),
week_errors as (select count(*) as week_errors from events where error_stack_trace is not null and event_time > current_timestamp - interval '7 days'),

day_events as (select count(*) as day_events from events where event_time > current_timestamp - interval '1 day'),
day_errors as (select count(*) as day_errors from events where error_stack_trace is not null and event_time > current_timestamp - interval '1 day')

select * from month_events, month_errors, week_events, week_errors, day_events, day_errors;
--end


--name: unique users by @date
--input: string date { pattern: datetime }
--connection: datagrok
select u.name, u.id from events e
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time)
group by u.name, u.id
--end


--name: events on @date
--input: string date { pattern: datetime }
--connection: datagrok
select u.login as user, u.id as user_id, t.name, t.source, e.event_time, e.description from events e
inner join event_types t on e.event_type_id = t.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time)
order by e.event_time desc
--end


--name: events on @date and @users
--input: string date { pattern: datetime }
--input: list users
--connection: datagrok
select u.login as user, u.id as user_id, t.name, t.source, e.event_time, e.description from events e
inner join event_types t on e.event_type_id = t.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time) and u.login = any(@users)
order by e.event_time desc
--end


--name: errors on @date
--input: string date { pattern: datetime }
--connection: datagrok
select friendly_name as name, event_time, error_message, error_stack_trace, source, run_number as count from events e
where @date(e.event_time) and e.error_stack_trace is not null
--end


--name: errors on @date and @users
--input: string date { pattern: datetime }
--input: list users
--connection: datagrok
select e.friendly_name as name, e.event_time, e.error_message, e.error_stack_trace, e.source, e.run_number as count from events e
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time) and e.error_stack_trace is not null and u.login = any(@users)
--end


--name: events summary on @date
--input: string date { pattern: datetime }
--connection: datagrok
select et.friendly_name, count(et.friendly_name) from events e
inner join event_types et on e.event_type_id = et.id
where @date(e.event_time)
group by et.friendly_name
order by count(et.friendly_name) desc
--end


--name: events summary on @date and @users
--input: string date { pattern: datetime }
--input: list users
--connection: datagrok
select et.friendly_name, count(et.friendly_name) from events e
inner join event_types et on e.event_type_id = et.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time) and u.login = any(@users)
group by et.friendly_name
order by count(et.friendly_name) desc
--end


--name: errors summary on @date
--input: string date { pattern: datetime }
--connection: datagrok
select concat(e.friendly_name, ': ', e.error_message), count(et.friendly_name) from events e
join event_types et on e.event_type_id = et.id
where @date(e.event_time) and e.error_stack_trace is not null
group by concat(e.friendly_name, ': ', e.error_message)
order by count(et.friendly_name) desc
--end


--name: errors summary on @date and @users
--input: string date { pattern: datetime }
--input: list users
--connection: datagrok
select concat(e.friendly_name, ': ', e.error_message), count(et.friendly_name) from events e
join event_types et on e.event_type_id = et.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time) and e.error_stack_trace is not null and u.login = any(@users)
group by concat(e.friendly_name, ': ', e.error_message)
order by count(et.friendly_name) desc
--end


--name: manual activity by @date
--input: string date { pattern: datetime }
select ta.date, ts.file_name, ta.result, ta.error_desc as error
from test_activity ta
join test_scenario ts on ts.id = ta.scenario_id
where @date(ta.date) and ts.type = 'manual'
order by ta.date desc
--end
