<!-- TITLE: Athena -->
<!-- SUBTITLE: -->

# Athena

Provides access to [Amazon Athena](https://docs.aws.amazon.com/athena/latest/ug/what-is.html) service
using SQL queries via JDBC driver . 

## Link CSVs from S3 with Athena

Simple steps to link CSVs with Athena and get results in the Datagrok:

* Upload CSVs to S3 bucket. Note that Athena looks into S3 folder not file, so if CSVs have different 
structure, they should be located in separate folders.

    For example: 
      * Bucket s3://athena-northwind/
        * orders/
            - orders.csv
        * products/
            - products.csv

* Create bucket or folder in existing bucket for Athena Output.

    For example: 
        s3://athena-northwind/results/

* Create table in [Athena console](https://console.aws.amazon.com/athena). UI builds SQL query for creating 
  table in Athena.

    Following example for Northwind "products.csv":
    
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
        * To skip header line add "skip.header.line.count" parameter

* Create [Data Connection](data-connection.md) in the Datagrok platforms.
    For example:
      - Name: northwind
      - Server: athena.us-east-2.amazonaws.com
      - Db: northwind
      - Port: 443
      - S3 Output Location: s3://athena-northwind/results/
      - Access Key: <key>
      - Secret Key: <secret>
      
    Notes:
      * Server has following format: "athena.<region>.amazonaws.com"
      * Do not forget "/" at "S3 Output Location" parameter end 
      
* Create [Data Query](data-query.md) under new connection:
    For example:
      
    ```sql
    SELECT * FROM northwind.products
    ```
 
    Notes:
       * Athena requires DB name before table name "northwind.products"
       
       
See also:

  * [Data Connection](../data-connection.md)
  * [Data Query](data-query.md)
  * [Amazon Athena](https://docs.aws.amazon.com/athena/latest/ug/what-is.html)
 