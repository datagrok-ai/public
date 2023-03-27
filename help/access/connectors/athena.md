---
title: "Athena"
---

Provides access to [Amazon Athena](https://docs.aws.amazon.com/athena/latest/ug/what-is.html)
service using SQL queries via a JDBC driver.

```json
{
   "region": "",
   "db": "",
   "S3OutputLocation": "",
   "VPCEndpoint":"",
   "S3OutputEncOption":""
}
```

## Supported Parameters

| Type                   | Value       | Description or Example     |
|------------------------|-------------|----------------------------|
| `num`, `int`, `double` | =           | =100                       |
|                        | >           | >1.02                      |
|                        | >=          | >=4.1                      |
|                        | <=          | <=100                      |
|                        | !=          | !=5                        |
|                        | in          | in (1, 3, 10.2)            |
|                        | min-max     | 1.5-10.0                   |
| `string`               | contains    | contains ea                |
|                        | starts with | starts with R              |
|                        | ends with   | ends with w                |
|                        | in          | in (ab, "c d", "e\\"f\\"") |
|                        | regex       | regex ^(.+)@(.+)$          |
| `datetime`             | anytime     |                            |
|                        | before      | before 1/1/2022            |
|                        | after       | after 1/1/2022             |
|                        | today       |                            |
|                        | this week   |                            |
|                        | this month  |                            |
|                        | this year   |                            |
|                        | last year   |                            |
|                        | min-max     |                            |
|                        | April 2021  |                            |
| `list<string>` (1)     |             |                            |

* (1) default parameters are not supported

## Supported output types

| Type                              | Supported              |
|-----------------------------------|------------------------|
| boolean                           | :white_check_mark:     |
| tinyint                           | :white_check_mark:     |
| smallint                          | :white_check_mark:     |
| int, integer                      | :white_check_mark:     |
| bigint                            | :white_check_mark:     |
| double                            | :white_check_mark:     |
| float                             | :white_check_mark:     |
| decimal                           | :white_check_mark:     |
| char, varchar, string             | :white_check_mark:     |
| date, timestamp                   | :white_check_mark:     |
| array<data_type>                  | :white_check_mark: (1) |
| map<primitive_type, data_type>    | :white_check_mark: (1) |
| struct<col_name : data_type  ...> | :white_check_mark: (1) |
| binary                            | not tested             |

* (1) supported as a string

## Supported features

* Schema browsing
* Build query
* Visual query
* Connection test

## Usage example: link csv files from s3 with athena

Simple steps to link csv files with Athena and get results in Datagrok:

1. Upload CSVs to an S3 bucket. Note that Athena looks into S3 folder, not file, so if CSVs have different structure,
   they should be located in separate folders. For example:

   ```
    Bucket s3://athena-northwind/
      orders/
        orders.csv
      products/
        products.csv
   ```

2. Create a bucket or folder in the existing bucket for Athena Output. For example:

   ```
   s3://athena-northwind/results/
   ```

3. Create table in [Athena console](https://console.aws.amazon.com/athena). UI builds SQL query for creating table in
   Athena. Following example for Northwind "products.csv":

   ```
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
   ```

   Tips:
   * To change CSV delimiter adjust "serialization.format" and "field.delim" parameters
   * To skip the header line add "skip.header.line.count" parameter

4. Create a [data connection](../data-connection.md) in the Datagrok platform. The <a href="#" id="parameters">
   parameters</a> may include: `"region"`, `"vpc endpoint"`, `"db"`
   , `"s3OutputLocation"`, `"s3OutputEncOption"`, `"accessKey"`, `"secretKey"`, or `"connString"`. For example:

   * Name: `northwind`
   * Region: `us-east-2`
   * Db: `northwind`
   * S3 Output Location: `s3://athena-northwind/results/`
   * Access Key: `<key>`
   * Secret Key: `<secret>`

     Notes:
      * VPC Endpoint is optional. If not specified, then canonical endpoint - "`athena.[Region].amazonaws.com:443`" will be used
      * Do not forget "/" at "S3 Output Location" parameter end

5. Create a [data query](../data-query.md) under the new connection. For example:

   ```
   SELECT * FROM products;
   ```

See also:

* [Data connection](../data-connection.md)
* [Data query](../data-query.md)
* [Amazon Athena](https://docs.aws.amazon.com/athena/latest/ug/what-is.html)
