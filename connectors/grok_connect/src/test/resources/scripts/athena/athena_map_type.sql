CREATE EXTERNAL TABLE `map_type`(
  `map_data` map<string,int>)
ROW FORMAT SERDE
  'org.apache.hadoop.hive.serde2.lazy.LazySimpleSerDe'
STORED AS INPUTFORMAT
  'org.apache.hadoop.mapred.TextInputFormat'
OUTPUTFORMAT
  'org.apache.hadoop.hive.ql.io.HiveIgnoreKeyTextOutputFormat'
LOCATION
  's3://datagrok-provider-test/tables/map-type'
TBLPROPERTIES (
  'transient_lastDdlTime'='1678120558')

INSERT INTO map_type VALUES(MAP(ARRAY['foo', 'bar', 'Datagrok'], ARRAY[1, 2, 2023]));
