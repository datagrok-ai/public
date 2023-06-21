CREATE TABLE date_types (date_data date, datetime_data datetime, datetime_data2 datetime2, time_data time,
datetimeoffset_data datetimeoffset, smalldatetime_data smalldatetime);

INSERT INTO date_types(date_data, datetime_data, datetime_data2, time_data, datetimeoffset_data, smalldatetime_data)
VALUES (CAST(N'1900-01-01' AS DATE), CAST(N'1900-01-01 00:00:00' AS DATETIME), CAST('2007-05-02T19:58:47.123' AS datetime2),
        CAST('00:00:00.00' AS TIME),
        CAST('1900-01-01 00:00:00 +02:00' AS datetimeoffset), CAST('1955-12-13 12:43:00' AS smalldatetime));