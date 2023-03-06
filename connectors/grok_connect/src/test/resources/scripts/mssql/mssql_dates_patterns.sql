SET DATEFIRST 1; -- set monday as first day of week (by default - sunday)

create table DATES_PATTERNS (
    date DATE
);

insert into DATES_PATTERNS (date) values (GETDATE());
insert into DATES_PATTERNS (date) values (DATEADD(day, -1, GETDATE())); -- yesterday

INSERT INTO DATES_PATTERNS(date) SELECT DATEADD(DAY, 7 - DATEPART(WEEKDAY, GETDATE()), GETDATE())
                                  WHERE NOT EXISTS (SELECT * FROM DATES_PATTERNS WHERE date = DATEADD(DAY, 7 - DATEPART(WEEKDAY, GETDATE()), GETDATE())); -- last day of week

insert into DATES_PATTERNS (date) values (DATEADD(day, -150, GETDATE()));
insert into DATES_PATTERNS (date) values ('2021-04-09');
