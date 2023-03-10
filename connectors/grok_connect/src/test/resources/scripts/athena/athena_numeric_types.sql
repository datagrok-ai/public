CREATE EXTERNAL TABLE `numeric_types`(
  `tinyint_data` tinyint,
  `smallint_data` smallint,
  `int_data` int,
  `bigint_data` bigint)
ROW FORMAT SERDE
  'org.apache.hadoop.hive.serde2.lazy.LazySimpleSerDe'
STORED AS INPUTFORMAT
  'org.apache.hadoop.mapred.TextInputFormat'
OUTPUTFORMAT
  'org.apache.hadoop.hive.ql.io.HiveIgnoreKeyTextOutputFormat'
LOCATION
  's3://datagrok-provider-test/tables/numeric-types'
TBLPROPERTIES (
  'transient_lastDdlTime'='1678111618')

INSERT INTO NUMERIC_TYPES (tinyint_data, smallint_data, int_data, bigint_data) VALUES (56, 1241, 2600000, 9223372036854775807);


INSERT INTO NUMERIC_TYPES (tinyint_data, smallint_data, int_data, bigint_data) VALUES (0, -1000, -2600000, -9223372036854775807);
