CREATE TABLE DATE_TYPES(date_type Date, date32_type Date32, datetime_type DateTime,
datetime64_type DateTime64(3, 'Europe/Kyiv'));

INSERT INTO DATE_TYPES(date_type, date32_type, datetime_type, datetime64_type)
VALUES (toDate('2016-06-15'), toDate('2016-06-15'), toDateTime('2016-06-15 23:00:00'),
        toDateTime('2019-01-01 00:00:00', 'UTC'));
