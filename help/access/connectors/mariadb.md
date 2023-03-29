---
title: "MariaDB"
---

Provides access to [MariaDB](https://mariadb.org/) database using SQL queries
via JDBC driver.

## Connection parameters

```json
{
  "server": "",
  "port": "",
  "db": "",
  "connString": ""
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

| Type                          | Supported              |
|-------------------------------|------------------------|
| INTEGER, SMALLINT             | :white_check_mark:     |
| DECIMAL, NUMERIC              | :white_check_mark:     |
| FLOAT, REAL, DOUBLE PRECISION | :white_check_mark:     |
| DATE, TIME                    | :white_check_mark:     |
| TIMESTAMP, YEAR               | :white_check_mark:     |
| CHAR, VARCHAR, TEXT           | :white_check_mark:     |
| JSON                          | :white_check_mark: (1) |
| GEOMETRY                      | limited support (2)    |
| BIT                           | limited support (3)    |
| BINARY, VARBINARY, BLOB       | not tested             |

* (1) supported as a string
* (2) you get unreadable representation, but in query you can cast such a types to text (e.g. ST_AsText function)
* (3) requires explicit cast (e.g. BIN function)

## Supported features

* Schema browsing
* Build query
* Visual query
* Connection test

See also:

* [Data connection](../data-connection.md)
* [MariaDB](https://mariadb.org/)
