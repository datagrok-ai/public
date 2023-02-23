CREATE TABLE date_patterns (
    dat DATE
);

INSERT INTO date_patterns(dat) VALUES (SYSDATE); --TODAY
INSERT INTO date_patterns(dat) VALUES (TRUNC(SYSDATE + 6, 'DAY') - 1); --THIS WEEK
INSERT INTO date_patterns(dat) VALUES (SYSDATE - 1); --YESTERDAY
INSERT INTO date_patterns(dat) VALUES (SYSDATE - 150);
INSERT INTO date_patterns(dat) VALUES (TO_DATE('2021-04-09', 'YYYY-MM-DD'));