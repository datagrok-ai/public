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
month_count as (select count(*) as "month" from users where joined > current_timestamp - interval '30 days'),
week_count as (select count(*) as "week" from users where joined > current_timestamp - interval '7 days'),
day_count as (select count(*) as "day" from users where joined > current_timestamp - interval '1 day')
select * from month_count, week_count, day_count;
--end


--name: new events last month week day
--connection: datagrok
with
month_events as (select count(*) as "month" from events where event_time > current_timestamp - interval '30 days'),
week_events as (select count(*) as "week" from events where event_time > current_timestamp - interval '7 days'),
day_events as (select count(*) as "day" from events where event_time > current_timestamp - interval '1 day')

select * from month_events, week_events, day_events;
--end

--name: new errors last month week day
--connection: datagrok
with
month_errors as (select count(*) as "month" from events where error_stack_trace is not null and event_time > current_timestamp - interval '30 days'),
week_errors as (select count(*) as "week" from events where error_stack_trace is not null and event_time > current_timestamp - interval '7 days'),
day_errors as (select count(*) as "day" from events where error_stack_trace is not null and event_time > current_timestamp - interval '1 day')

select * from month_errors, week_errors, day_errors;
--end


--name: unique users per day on @date and @users and @events
--input: string date { pattern: datetime }
--input: list users
--input: list events
--input: bool isExactly
--connection: datagrok
select t.date::date, count(t.id) as user_count from (
	select distinct on (date(e.event_time), u.id) u.id, e.event_time as date
	from events e
	inner join users_sessions s on e.session_id = s.id
	inner join users u on u.id = s.user_id
	where @date(e.event_time)
	and (u.login = any(@users) or @users = ARRAY['all'])
    and
    case @isExactly
    when true then (e.friendly_name = any(@events) or @events = ARRAY['all'])
    else (e.friendly_name like any (@events) or @events = ARRAY['all'])
    end
) t
group by t.date::date;
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


----name: events on @date
----input: string date { pattern: datetime }
----connection: datagrok
--select u.login as user, u.id as user_id, t.name, t.source, e.id as event_id, e.event_time, e.description from events e
--inner join event_types t on e.event_type_id = t.id
--inner join users_sessions s on e.session_id = s.id
--inner join users u on u.id = s.user_id
--where @date(e.event_time)
--order by e.event_time desc
----end


--name: events on @date and @users and @events
--input: string date { pattern: datetime }
--input: list users
--input: list events
--input: bool isExactly
--connection: datagrok
select u.login as user, u.id as user_id, t.name, t.source, e.id as event_id, e.event_time, e.description from events e
inner join event_types t on e.event_type_id = t.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time)
and (u.login = any(@users) or @users = ARRAY['all'])
and
case @isExactly
when true then (e.friendly_name = any(@events) or @events = ARRAY['all'])
else (e.friendly_name like any (@events) or @events = ARRAY['all'])
end
order by e.event_time desc
--end


----name: errors on @date
----input: string date { pattern: datetime }
----connection: datagrok
--select friendly_name as name, event_time, error_message, error_stack_trace, source, run_number as count, e.id as event_id from events e
--where @date(e.event_time) and e.error_stack_trace is not null
----end


--name: errors on @date and @users and @events
--input: string date { pattern: datetime }
--input: list users
--input: list events
--input: bool isExactly
--connection: datagrok
select e.friendly_name as name, e.event_time, e.error_message, e.error_stack_trace, e.source, e.run_number as count, e.id as event_id from events e
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where e.error_stack_trace is not null
and @date(e.event_time)
and (u.login = any(@users) or @users = ARRAY['all'])
and
case @isExactly
when true then (e.friendly_name = any(@events) or @events = ARRAY['all'])
else (e.friendly_name like any (@events) or @events = ARRAY['all'])
end
--end


--name: events summary on @date
--input: string date { pattern: datetime }
--connection: datagrok
select et.friendly_name, et.id as event_type_id, count(1) as cnt from events e
inner join event_types et on e.event_type_id = et.id
where @date(e.event_time)
group by et.friendly_name, et.id
order by cnt desc
limit 25
--end


--name: events summary on @date and @users
--input: string date { pattern: datetime }
--input: list users
--connection: datagrok
select et.friendly_name, et.id as event_type_id, count(1) as cnt from events e
inner join event_types et on e.event_type_id = et.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time) and u.login = any(@users)
group by et.friendly_name, et.id
order by cnt desc
limit 25
--end


--name: errors summary on @date
--input: string date { pattern: datetime }
--connection: datagrok
select et.friendly_name, et.id as event_type_id, count(1) as cnt from events e
join event_types et on e.event_type_id = et.id
where @date(e.event_time) and e.error_stack_trace is not null
group by et.friendly_name, et.id
order by cnt desc
limit 25
--end


--name: errors summary on @date and @users
--input: string date { pattern: datetime }
--input: list users
--connection: datagrok
select et.friendly_name, et.id as event_type_id, count(1) as cnt from events e
join event_types et on e.event_type_id = et.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time) and e.error_stack_trace is not null and u.login = any(@users)
group by et.friendly_name, et.id
order by cnt desc
limit 25
--end


--name: manual activity by @date
--input: string date { pattern: datetime }
--connection: TestTrack
select ta.date, ts.file_name, ta.result, ta.error_desc as error
from test_activity ta
join test_scenario ts on ts.id = ta.scenario_id
where @date(ta.date) and ts.type = 'manual'
order by ta.date desc
--end
