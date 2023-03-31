---
title: "DB2"
---

Provides access to [IBM Db2](https://www.ibm.com/analytics/db2) service using
SQL queries via JDBC driver.

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

| Type                                        | Supported              |
|---------------------------------------------|------------------------|
| INTEGER, SMALLINT, BIGINT                   | :white_check_mark:     |
| BOOLEAN                                     | :white_check_mark:     |
| DECIMAL, DECFLOAT                           | :white_check_mark:     |
| REAL, DOUBLE                                | :white_check_mark:     |
| DATE, TIME, TIMESTAMP                       | :white_check_mark:     |
| CHARACTER, VARCHAR, GRAPHIC,<br/>VARGRAPHIC | :white_check_mark:     |
| CLOB, DBCLOB                                | :white_check_mark:     |
| XML                                         | :white_check_mark: (1) |
| BINARY, VARBINARY, BLOB                     | not tested             |

* (1) supported as a string

## Supported features

* Schema browsing
* Build query
* Visual query
* Connection test

See also:

* [Data connection](../data-connection.md)
* [IBM Db2](https://www.ibm.com/analytics/db2)
* [IBM Db2 Family](https://en.wikipedia.org/wiki/IBM_Db2_Family)
