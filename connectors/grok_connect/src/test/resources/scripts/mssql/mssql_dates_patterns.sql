create table DATES_PATTERNS (
    "date" DATE
);

insert into DATES_PATTERNS ("date")values (GETDATE());
insert into DATES_PATTERNS ("date")values (DATEADD(day, -1, GETDATE())); -- yesterday
insert into DATES_PATTERNS ("date")values (DATEADD(DAY, 8 - DATEPART(WEEKDAY, GETDATE()), GETDATE())); -- last day of this week
insert into DATES_PATTERNS ("date")values (DATEADD(day, -150, GETDATE()));
insert into DATES_PATTERNS ("date")values ('2021-04-09');