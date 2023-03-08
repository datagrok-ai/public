CREATE EXTERNAL TABLE `float_types`(
  `double_type` double,
  `float_type` float,
  `decimal_type` decimal(11,2))
ROW FORMAT SERDE
  'org.apache.hadoop.hive.serde2.lazy.LazySimpleSerDe'
STORED AS INPUTFORMAT
  'org.apache.hadoop.mapred.TextInputFormat'
OUTPUTFORMAT
  'org.apache.hadoop.hive.ql.io.HiveIgnoreKeyTextOutputFormat'
LOCATION
  's3://datagrok-provider-test/tables/float-type'
TBLPROPERTIES (
  'transient_lastDdlTime'='1678118970')

INSERT INTO float_types (double_type, float_type, decimal_type) VALUES (1.79769313486231570e+308,
                                                                        3.40282346638528860e+38, 0.12);

INSERT INTO float_types (double_type, float_type, decimal_type) VALUES (4.94065645841246544e-324,
                                                                        1.40129846432481707e-45, 0.50);
