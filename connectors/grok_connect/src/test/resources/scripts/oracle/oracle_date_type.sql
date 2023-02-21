CREATE TABLE dates_type (
    date DATE,
    stamp TIMESTAMP,
    zoned_stamp TIMESTAMP WITH TIME ZONE,
    interval1 INTERVAL YEAR TO MONTH,
    interval2 INTERVAL DAY TO SECOND
);

INSERT INTO dates_type("date", stamp, zoned_stamp, interval1, interval2)
VALUES ('01-01-2023', '03-AUG-17 11:20:30.45 AM',
        TO_TIMESTAMP_TZ ('21-FEB-2009 18:00:00 -5:00', 'DD-MON-YYYY HH24:MI:SS TZH:TZM'),
        INTERVAL '10-2' YEAR TO MONTH, INTERVAL '4 5:12:10.222' DAY TO SECOND);