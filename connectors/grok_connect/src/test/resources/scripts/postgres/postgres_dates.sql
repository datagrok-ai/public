CREATE TABLE DATES(date DATE, time TIME, stamp TIMESTAMP, interval INTERVAL);

INSERT INTO DATES("date", "time", stamp, "interval") VALUES ('1999-01-08', '04:05:06.789', '1999-01-08 04:05:06',
                                                             '1 year 5 months 5 days');
INSERT INTO DATES("date", "time", stamp, "interval") VALUES ('January 8, 1999', '04:05:06', '2004-10-19 10:23:54',
                                                             '1 day');
INSERT INTO DATES("date", "time", stamp, "interval") VALUES ('1/8/1999', '040506', '2004-10-19 10:23:54+02',
                                                             '1 second');
INSERT INTO DATES("date", "time", stamp, "interval") VALUES ('1/18/1999', '04:05 AM',
                                                             '2004-10-19 10:23:54',
                                                             '5 years 4 months 3 days 2 hours 1 minute 1 second');
