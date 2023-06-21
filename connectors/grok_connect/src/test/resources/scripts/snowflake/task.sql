create or replace task TEST.PUBLIC.UPDATE_DATES_PATTERNS
	warehouse=COMPUTE_WH
	schedule='USING CRON 0 2 * * * Europe/Kiev'
	as execute immediate
$$
BEGIN
DELETE FROM dates_patterns;
insert into DATES_PATTERNS (date) values (CURRENT_DATE);
insert into DATES_PATTERNS (date) values (CURRENT_DATE - 1); -- yesterday
insert into DATES_PATTERNS (date) select (LAST_DAY(CURRENT_DATE, 'WEEK'))
    where not exists (select * from DATES_PATTERNS where date = LAST_DAY(CURRENT_DATE, 'WEEK')); -- last day of this week
insert into DATES_PATTERNS (date) values (CURRENT_DATE - 150);
insert into DATES_PATTERNS (date) values ('2021-04-09');
END;
$$;
