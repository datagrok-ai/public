CREATE EXTERNAL TABLE IF NOT EXISTS northwind.products (
  `productid` int,
  `productname` string,
  `supplierid` int,
  `categoryid` int,
  `quantityperunit` string,
  `unitprice` double,
  `unitsinstock` int,
  `unitsonorder` int,
  `reorderlevel` int,
  `discontinued` int
)
ROW FORMAT SERDE 'org.apache.hadoop.hive.serde2.lazy.LazySimpleSerDe'
WITH SERDEPROPERTIES (
  'serialization.format' = ',',
  'field.delim' = ','
) LOCATION 's3://athena-northwind/products/'
TBLPROPERTIES (
  'has_encrypted_data' = 'false',
  'skip.header.line.count' = '1'
);


CREATE EXTERNAL TABLE IF NOT EXISTS northwind.orders (
  `orderid` int,
  `customerid` string,
  `employeeid` int,
  `orderdate` timestamp,
  `requireddate` timestamp,
  `shippeddate` timestamp,
  `shipvia` int,
  `freight` double,
  `shipname` string,
  `shipaddress` string,
  `shipcity` string,
  `shipregion` string,
  `shippostalcode` int,
  `shipcountry` string
)
ROW FORMAT SERDE 'org.apache.hadoop.hive.serde2.lazy.LazySimpleSerDe'
WITH SERDEPROPERTIES (
  'serialization.format' = ',',
  'field.delim' = ','
) LOCATION 's3://athena-northwind/orders/'
TBLPROPERTIES (
  'has_encrypted_data' = 'false',
  'skip.header.line.count' = '1'
);
