DELETE FROM dates_patterns;

insert into dates_patterns (date) values (CURDATE());
insert into dates_patterns (date) values (DATE_SUB(CURDATE(), INTERVAL 1 DAY)); -- yesterday

INSERT INTO dates_patterns(date) SELECT DATE(NOW() + INTERVAL (6 - WEEKDAY(NOW())) DAY) AS date
FROM dates_patterns WHERE (date = DATE(NOW() + INTERVAL (6 - WEEKDAY(NOW())) DAY)) HAVING COUNT(*) = 0; -- last day of week

insert into dates_patterns (date) values (DATE_SUB(CURDATE(), INTERVAL 150 DAY));
insert into dates_patterns (date) values ('2021-04-09');
