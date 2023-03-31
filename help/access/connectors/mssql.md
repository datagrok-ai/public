---
title: "MS SQL"
---

Provides access to [Microsoft SQL](https://www.microsoft.com/en-us/sql-server) database using SQL queries via JDBC
driver.

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

| Type                                                                | Supported              |
|---------------------------------------------------------------------|------------------------|
| bigint                                                              | :white_check_mark:     |
| bit                                                                 | :white_check_mark:     |
| decimal                                                             | :white_check_mark:     |
| int                                                                 | :white_check_mark:     |
| smallint                                                            | :white_check_mark:     |
| tinyint                                                             | :white_check_mark:     |
| numeric                                                             | :white_check_mark:     |
| money                                                               | :white_check_mark:     |
| decimal                                                             | :white_check_mark:     |
| float                                                               | :white_check_mark:     |
| real                                                                | :white_check_mark:     |
| date, datetime, datetime2, time, <br/>datetimeoffset, smalldatetime | :white_check_mark:     |
| xml                                                                 | :white_check_mark: (1) |
| geography, geometry                                                 | :white_check_mark: (1) |
| xml                                                                 | :white_check_mark: (1) |
| binary, varbinary                                                   | limited support    (2) |

* (1) supported as a string
* (2) you get unreadable representation, but in query you can cast such a types to varchar, hex

## Supported features

* Schema browsing
* Build query
* Visual query
* Connection test

See also:

* [Data connection](../data-connection.md)
* [Microsoft SQL](https://www.microsoft.com/en-us/sql-server)
* [Microsoft SQL Wiki](https://en.wikipedia.org/wiki/Microsoft_SQL_Server)
