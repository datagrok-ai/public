CREATE TABLE DATE_TYPES(date_type date, time_type time, timetz_type timetz, timestamp_type timestamp,
timestamptz_type timestamptz);

INSERT INTO DATE_TYPES(date_type, time_type, timetz_type, timestamp_type, timestamptz_type)
VALUES ('2008-06-01', '00:00:00', '04:05:06.789', '2008-06-01 09:59:59', '2001-02-16 18:38:40');
