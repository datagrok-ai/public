create table DATES_PATTERNS (
    date DATE
);

insert into DATES_PATTERNS (date) values (current_date);
insert into DATES_PATTERNS (date) values (current_date - 1); -- yesterday
insert into DATES_PATTERNS (date) values ((date_trunc('week', current_date::timestamp) + '6 days'::interval)::date); -- last day of this week
insert into DATES_PATTERNS (date) values (current_date - 150);
insert into DATES_PATTERNS (date) values ('2021-04-09');
