--name: NewUsersEventsErrors
--connection: System:DatagrokAdmin
with
month_count as (select count(*) as "month" from users where joined > current_timestamp - interval '30 days'),
week_count as (select count(*) as "week" from users where joined > current_timestamp - interval '7 days'),
day_count as (select count(*) as "day" from users where joined > current_timestamp - interval '1 day'),
month_events as (select count(*) as "month" from events where event_time > current_timestamp - interval '30 days'),
week_events as (select count(*) as "week" from events where event_time > current_timestamp - interval '7 days'),
day_events as (select count(*) as "day" from events where event_time > current_timestamp - interval '1 day'),
month_errors as (select count(*) as "month" from events where error_stack_trace is not null and event_time > current_timestamp - interval '30 days'),
week_errors as (select count(*) as "week" from events where error_stack_trace is not null and event_time > current_timestamp - interval '7 days'),
day_errors as (select count(*) as "day" from events where error_stack_trace is not null and event_time > current_timestamp - interval '1 day')
select 'users_count' as counter_type, * from month_count, week_count, day_count union
select 'events_count' as counter_type, * from month_events, week_events, day_events union
select 'errors_count' as counter_type, * from month_errors, week_errors, day_errors
--end


--name: UniqueUsersPerDayOnDateAndUsersAndEvents
--input: string date { pattern: datetime }
--input: list users
--input: list events
--input: bool isExactly
--connection: System:DatagrokAdmin
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


--name: UniqueUsersByDate
--input: string date { pattern: datetime }
--connection: System:DatagrokAdmin
select u.name, u.id from events e
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time)
group by u.name, u.id
--end

--name: EventsOnDateAndUsersAndEvents
--input: string date { pattern: datetime }
--input: list users
--input: list events
--input: bool isExactly
--connection: System:DatagrokAdmin
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


--name: ErrorsOnDateAndUsersAndEvents
--input: string date { pattern: datetime }
--input: list users
--input: list events
--input: bool isExactly
--connection: System:DatagrokAdmin
select et.friendly_name, count(e.id) as count, et.error_source, et.error_stack_trace, et.id from event_types et
join events e on e.event_type_id = et.id
join users_sessions s on e.session_id = s.id
join users u on u.id = s.user_id
where et.source = 'error'
and @date(e.event_time)
and (u.login = any(@users) or @users = ARRAY['all'])
and
case @isExactly
when true then (e.friendly_name = any(@events) or @events = ARRAY['all'])
else (e.friendly_name like any (@events) or @events = ARRAY['all'])
end
group by et.id, et.friendly_name, et.error_stack_trace, et.error_source
order by count desc
--end


--name: EventsSummaryOnDate
--input: string date { pattern: datetime }
--connection: System:DatagrokAdmin
select et.friendly_name, et.id as event_type_id, count(1) as cnt from events e
inner join event_types et on e.event_type_id = et.id
where @date(e.event_time)
group by et.friendly_name, et.id
order by cnt desc
limit 25
--end


--name: EventsSummaryOnDateAndUsers
--input: string date { pattern: datetime }
--input: list users
--connection: System:DatagrokAdmin
select et.friendly_name, et.id as event_type_id, count(1) as cnt from events e
inner join event_types et on e.event_type_id = et.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time) and u.login = any(@users)
group by et.friendly_name, et.id
order by cnt desc
limit 25
--end


--name: ErrorsSummaryOnDate
--input: string date { pattern: datetime }
--connection: System:DatagrokAdmin
select et.friendly_name, et.id as event_type_id, count(1) as cnt from events e
join event_types et on e.event_type_id = et.id
where @date(e.event_time) and e.error_stack_trace is not null
group by et.friendly_name, et.id
order by cnt desc
limit 25
--end


--name: ErrorsSummaryOnDateAndUsers
--input: string date { pattern: datetime }
--input: list users
--connection: System:DatagrokAdmin
select et.friendly_name, et.id as event_type_id, count(1) as cnt from events e
join event_types et on e.event_type_id = et.id
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where @date(e.event_time) and e.error_stack_trace is not null and u.login = any(@users)
group by et.friendly_name, et.id
order by cnt desc
limit 25
--end


--name: errors for @user
--input: string date { pattern: datetime }
--input: string user
--connection: System:DatagrokAdmin
select e.friendly_name as name, e.event_time, e.error_message, e.error_stack_trace, e.source, e.run_number as count, e.id as event_id from events e
inner join users_sessions s on e.session_id = s.id
inner join users u on u.id = s.user_id
where e.error_stack_trace is not null
and @date(e.event_time)
and (u.login = @user)
--end