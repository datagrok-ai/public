CREATE EXTERNAL TABLE `array_type_int`(
  `id` int,
  `array_data` array<int>)
ROW FORMAT SERDE
  'org.apache.hadoop.hive.serde2.lazy.LazySimpleSerDe'
STORED AS INPUTFORMAT
  'org.apache.hadoop.mapred.TextInputFormat'
OUTPUTFORMAT
  'org.apache.hadoop.hive.ql.io.HiveIgnoreKeyTextOutputFormat'
LOCATION
  's3://datagrok-provider-test/tables/array-type-int'
TBLPROPERTIES (
  'transient_lastDdlTime'='1678120400')

INSERT INTO array_type_int values(1, ARRAY[1,2,3]);
INSERT INTO array_type_int values(2, ARRAY[1,NULL,3]);

CREATE EXTERNAL TABLE `array_type_string`(
  `id` int,
  `array_data` array<string>)
ROW FORMAT SERDE
  'org.apache.hadoop.hive.serde2.lazy.LazySimpleSerDe'
STORED AS INPUTFORMAT
  'org.apache.hadoop.mapred.TextInputFormat'
OUTPUTFORMAT
  'org.apache.hadoop.hive.ql.io.HiveIgnoreKeyTextOutputFormat'
LOCATION
  's3://datagrok-provider-test/tables/array-type-str'
TBLPROPERTIES (
  'transient_lastDdlTime'='1678120334')

INSERT INTO array_type_string values(1, ARRAY['Hello', 'World', 'Datagrok']);
