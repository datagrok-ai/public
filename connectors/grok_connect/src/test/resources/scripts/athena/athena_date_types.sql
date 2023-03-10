CREATE EXTERNAL TABLE `date_types`(
    `id` int,
  `date_type` date,
  `timestamp_type` timestamp)
ROW FORMAT SERDE
  'org.apache.hadoop.hive.serde2.lazy.LazySimpleSerDe'
STORED AS INPUTFORMAT
  'org.apache.hadoop.mapred.TextInputFormat'
OUTPUTFORMAT
  'org.apache.hadoop.hive.ql.io.HiveIgnoreKeyTextOutputFormat'
LOCATION
  's3://datagrok-provider-test/tables/date-types'
TBLPROPERTIES (
  'transient_lastDdlTime'='1678119775')


INSERT INTO date_types VALUES (1, DATE'1996-08-26', TIMESTAMP'2018-04-01 00:00:00');
INSERT INTO date_types VALUES (2, DATE'2023-08-26', TIMESTAMP'2023-04-05 12:00:00.123');
