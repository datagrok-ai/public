CREATE EXTERNAL TABLE `character_types`(
  `char_data` char(20),
  `varchar_data` varchar(20),
  `string_data` string)
ROW FORMAT SERDE
  'org.apache.hadoop.hive.serde2.lazy.LazySimpleSerDe'
STORED AS INPUTFORMAT
  'org.apache.hadoop.mapred.TextInputFormat'
OUTPUTFORMAT
  'org.apache.hadoop.hive.ql.io.HiveIgnoreKeyTextOutputFormat'
LOCATION
  's3://datagrok-provider-test/tables/character-types'
TBLPROPERTIES (
  'transient_lastDdlTime'='1678119322')

INSERT INTO character_types VALUES ('Datagrok', 'Hello, World', 'Hello, Datagrok!');
