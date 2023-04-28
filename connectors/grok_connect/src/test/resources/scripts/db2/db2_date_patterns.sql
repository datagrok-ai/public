CREATE TABLE dates_patterns (date DATE);

INSERT INTO dates_patterns(date) VALUES (CURRENT DATE);
INSERT INTO dates_patterns(date) VALUES (CURRENT DATE - 1 DAY);
INSERT INTO dates_patterns(date) SELECT (CURRENT DATE + (7 - DAYOFWEEK_ISO(CURRENT DATE)) DAYS) FROM sysibm.sysdummy1 WHERE NOT EXISTS (SELECT * FROM dates_patterns WHERE date = CURRENT DATE + (7 - DAYOFWEEK_ISO(CURRENT DATE)) DAYS);
INSERT INTO dates_patterns(date) VALUES (CURRENT DATE - 150 DAYS);
INSERT INTO dates_patterns(date) VALUES ('2021-04-09');
