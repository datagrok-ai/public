CREATE EXTERNAL TABLE `dates_patterns`(
  `date_data` date)
ROW FORMAT SERDE
  'org.apache.hadoop.hive.serde2.lazy.LazySimpleSerDe'
STORED AS INPUTFORMAT
  'org.apache.hadoop.mapred.TextInputFormat'
OUTPUTFORMAT
  'org.apache.hadoop.hive.ql.io.HiveIgnoreKeyTextOutputFormat'
LOCATION
  's3://datagrok-provider-test/tables/dates-patterns'
TBLPROPERTIES (
  'transient_lastDdlTime'='1678121587')

INSERT INTO dates_patterns VALUES (current_date); --today
INSERT INTO dates_patterns VALUES (date_add('day', -1, current_date)); --yesterday
INSERT INTO dates_patterns SELECT date_add('day', 6, date_trunc('week', current_date)) WHERE
NOT EXISTS (SELECT * FROM dates_patterns WHERE date_data = date_add('day', 6, date_trunc('week', current_date)));
INSERT INTO dates_patterns VALUES (date_add('day', -150, current_date));
INSERT INTO dates_patterns VALUES (DATE'2021-04-09');
