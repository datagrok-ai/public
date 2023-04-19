---
title: "Teradata"
---

Provides access to
[Teradata](https://www.teradata.ru/Products/Software/Database) database using
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

| Type                                                                         | Supported              |
|------------------------------------------------------------------------------|------------------------|
| BYTEINT, SMALLINT, INTEGER                                                   | :white_check_mark:     |
| BIGINT                                                                       | :white_check_mark:     |
| FLOAT, REAL, DECIMAL                                                         | :white_check_mark:     |
| CHAR, VARCHAR, CLOB                                                          | :white_check_mark:     |
| DATE, TIME, TIMESTAMP, TIME <br/>WITH TIMEZONE, <br/>TIMESTAMP WITH TIMEZONE | :white_check_mark:     |
| JSON                                                                         | :white_check_mark: (1) |
| XML                                                                          | :white_check_mark: (1) |
| SPATIAL                                                                      | :white_check_mark: (1) |
| ARRAY                                                                        | :white_check_mark: (1) |
| BYTE, VARBYTE, BLOB                                                          | limited support    (2) |

* (1) supported as a string
* (2) you get unreadable representation, but in query you can cast such a types to varchar, hex

## Supported features

* Schema browsing
* Build query
* Visual query
* Connection test

See also:

* [Data connection](../data-connection.md)
* [Teradata](https://www.teradata.ru/Products/Software/Database)
