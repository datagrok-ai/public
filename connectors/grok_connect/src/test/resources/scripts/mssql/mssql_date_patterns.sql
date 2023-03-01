create table DATES_PATTERNS (
    date_data DATE
);

insert into DATES_PATTERNS (date_data) values (CONVERT(date, GETDATE());
insert into DATES_PATTERNS (date_data) values (DATEADD(day, -1, CONVERT(date, GETDATE())); -- yesterday
insert into DATES_PATTERNS (date_data) values (DATEADD(DAY, 7 - DATEPART(WEEKDAY, GETDATE()), CAST(GETDATE() AS DATE))); -- last day of this week
insert into DATES_PATTERNS (date_data) values (DATEADD(day, -150, CONVERT(date, GETDATE()));
insert into DATES_PATTERNS (date_data) values ('2021-04-09');