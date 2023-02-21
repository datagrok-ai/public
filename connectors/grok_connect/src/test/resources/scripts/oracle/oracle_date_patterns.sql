CREATE TABLE date_patterns (
    date DATE
);

INSERT INTO date_patterns(date) VALUES (SYSDATE); --TODAY
INSERT INTO date_patterns(date) VALUES (TRUNC(SYSDATE + 6, 'DAY') - 1); --THIS WEEK
INSERT INTO date_patterns(date) VALUES (SYSDATE - 1); --YESTERDAY
INSERT INTO date_patterns(date) VALUES (TRUNC(SYSDATE +150, 'YEAR') - 1); --THIS YEAR
INSERT INTO date_patterns(date) VALUES ('2021-04-09');