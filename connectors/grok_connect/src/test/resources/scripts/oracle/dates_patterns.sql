DELETE FROM dates_patterns;

INSERT INTO dates_patterns(dat) VALUES (SYSDATE); --TODAY
INSERT INTO dates_patterns(dat) VALUES (SYSDATE - 1); --YESTERDAY
INSERT INTO dates_patterns(dat) VALUES (TRUNC(SYSDATE + 6, 'DAY')); --THIS WEEK
INSERT INTO dates_patterns(dat) VALUES (SYSDATE - 150);
INSERT INTO dates_patterns(dat) VALUES (TO_DATE('2021-04-09', 'YYYY-MM-DD'));