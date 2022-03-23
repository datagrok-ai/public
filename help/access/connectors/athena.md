<!-- TITLE: Athena -->
<!-- SUBTITLE: -->

# Athena

Provides access to [Amazon Athena](https://docs.aws.amazon.com/athena/latest/ug/what-is.html)
service using SQL queries via a JDBC driver.

## Link csvs from s3 with athena

Simple steps to link CSVs with Athena and get results in Datagrok:

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

   ```sql
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
   parameters</a> may include: `"server"`, `"db"`, `"port"`
   , `"s3OutputLocation"`, `"s3OutputEncOption"`, `"accessKey"`, `"secretKey"`, or `"connString"`. For example:

    * Name: `northwind`
    * Server: `athena.us-east-2.amazonaws.com`
    * Db: `northwind`
    * Port: `443`
    * S3 Output Location: `s3://athena-northwind/results/`
    * Access Key: `<key>`
    * Secret Key: `<secret>`

      Notes:
        * Server has the following format: "`athena.<region>.amazonaws.com`"
        * Do not forget "/" at "S3 Output Location" parameter end

5. Create a [data query](../data-query.md) under the new connection. For example:

   ```sql
   SELECT * FROM northwind.products
   ```

   Notes:
    * Athena requires DB name before table name "northwind.products"

See also:

* [Data connection](../data-connection.md)
* [Data query](../data-query.md)
* [Amazon Athena](https://docs.aws.amazon.com/athena/latest/ug/what-is.html)
