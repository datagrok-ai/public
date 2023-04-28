CREATE TABLE null_safety(id bigint, varchar_data varchar(max), char_data char(20), binary_data binary(8),
varbinary_data varbinary(6), text_data text, nchar_data nchar(16), nvarchar_data nvarchar(max), ntext_data ntext,
date_data date, datetime_data datetime, datetime_data2 datetime2, time_data time, datetimeoffset_data datetimeoffset,
smalldatetime_data smalldatetime, money_type money, small_money smallmoney, int_data int, smallint_data smallint, tinyint_data tinyint, bit_data bit,
    decimal_data decimal(5, 2), numeric_data numeric(10, 5), GeomCol1 geometry, xml_data xml
);

INSERT INTO null_safety(id, varchar_data, char_data, binary_data, varbinary_data, text_data, nchar_data,
                        nvarchar_data, ntext_data, date_data, datetime_data, datetime_data2, time_data,
                        datetimeoffset_data, smalldatetime_data, money_type, small_money, int_data, smallint_data,
                        tinyint_data, bit_data, decimal_data, numeric_data, GeomCol1, xml_data)
VALUES (NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
        NULL, NULL, NULL, NULL, NULL, NULL, NULL);
