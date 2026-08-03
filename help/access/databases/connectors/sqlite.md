---
title: "SQLite"
description: Connect Datagrok to a SQLite database file stored on a file share and query it with SQL via a JDBC driver.
keywords:
  - sqlite file
  - embedded database
  - db file
  - jdbc driver
  - connect to a database
  - file share
---

Provides access to [SQLite](https://www.sqlite.org/index.html) DB file using SQL
queries via a JDBC driver. The file should be stored on a Datagrok file share for
SQLite DBs.

## Connection parameters

```json
{
  "db": "",
  "connString": ""
}
```

See also:

* [Data connection](../../access.md#data-connection)
* [SQLite](https://www.sqlite.org/index.html)
