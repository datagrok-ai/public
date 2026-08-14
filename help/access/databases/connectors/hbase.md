---
title: "HBase"
description: Connect Datagrok to an Apache HBase NoSQL database and query it using SQL over JDBC.
keywords:
  - apache hbase
  - hbase nosql database
  - connect to hbase
  - hbase jdbc driver
  - sql queries over nosql
  - wide-column store
---

Provides access to [Apache HBase](https://hbase.apache.org/) NoSQL database
using SQL queries via a JDBC driver.

:::warning
The HBase connector is currently unavailable: the platform does not ship its JDBC (Phoenix
query server) driver, so the connector is not advertised. Contact Datagrok if you need it.
:::

## Connection parameters

```json
{
  "server": "",
  "port": "",
  "db": "",
  "connString": ""
}
```

See also:

* [Data connection](../../access.md#data-connection)
* [Apache HBase](https://hbase.apache.org/)
* [Apache HBase Wiki](https://en.wikipedia.org/wiki/Apache_HBase)
