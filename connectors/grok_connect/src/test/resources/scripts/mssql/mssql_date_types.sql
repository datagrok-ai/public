CREATE TABLE date_types (date_data date, datetime_data datetime, datetime_data2 datetime2, time_data time,
datetimeoffset_data datetimeoffset, smalldatetime_data smalldatetime);

INSERT INTO date_types(date_data, datetime_data, datetime_data2, time_data, datetimeoffset_data, smalldatetime_data)
VALUES ('1900-01-01', '1900-01-01 00:00:00', '2007-05-02T19:58:47.1234567', '23:59:59.9999999',
        '1900-01-01 00:00:00 00:00', '1955-12-13 12:43:00');