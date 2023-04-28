CREATE TABLE dates_patterns (
    dat DATE
);

INSERT INTO dates_patterns(dat) VALUES (SYSDATE); --TODAY
INSERT INTO dates_patterns(dat) VALUES (SYSDATE - 1); --YESTERDAY
INSERT INTO dates_patterns(dat) SELECT TRUNC(SYSDATE + 6, 'DAY') FROM dual
WHERE NOT EXISTS (SELECT * FROM dates_patterns WHERE to_date(dat) = to_date(TRUNC(SYSDATE + 6, 'DAY')));
INSERT INTO dates_patterns(dat) VALUES (SYSDATE - 150);
INSERT INTO dates_patterns(dat) VALUES (TO_DATE('2021-04-09', 'YYYY-MM-DD'));
