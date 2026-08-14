---
title: "PI"
description: Connect Datagrok to an AVEVA PI System process historian and query it with SQL via a JDBC driver.
keywords:
  - aveva pi system
  - process historian
  - pi database
  - jdbc driver
  - connect to a database
---

Provides access to the
[PI](https://www.aveva.com/en/products/aveva-pi-system/) database using SQL queries via a JDBC driver.

## Connection parameters

````json
{
   "parameters": {
    "accessServer": "",
    "port": port,
    "db": "db-name",
    "cacheSchema": false,
    "cacheResults": true,
    "ssl": false,
    "connString": ""
  }
}
````

See also:

* [Data connection](../../access.md#data-connection)
* [PI](https://www.aveva.com/en/products/aveva-pi-system/)
