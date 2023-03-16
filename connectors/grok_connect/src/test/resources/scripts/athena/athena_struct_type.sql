CREATE EXTERNAL TABLE `struct_type`(
  `struct_data` struct<name:varchar(10),age:int>)
ROW FORMAT SERDE
  'org.apache.hadoop.hive.serde2.lazy.LazySimpleSerDe'
STORED AS INPUTFORMAT
  'org.apache.hadoop.mapred.TextInputFormat'
OUTPUTFORMAT
  'org.apache.hadoop.hive.ql.io.HiveIgnoreKeyTextOutputFormat'
LOCATION
  's3://datagrok-provider-test/tables/struct-type'
TBLPROPERTIES (
  'transient_lastDdlTime'='1678120995')

INSERT INTO struct_type SELECT CAST(ROW('Bob', 38) AS ROW(name VARCHAR(10), age INTEGER));
