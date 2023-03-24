---
title: "Oracle"
---

Provides access to [Oracle Database](https://www.oracle.com/database/) database
using SQL queries via the JDBC driver.

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

| Type            | Supported              |
|-----------------|------------------------|
| NUMBER          | :white_check_mark:     |
| FLOAT           | :white_check_mark:     |
| VARCHAR2        | :white_check_mark:     |
| NVARCHAR2       | :white_check_mark:     |
| CHAR            | :white_check_mark:     |
| NCHAR           | :white_check_mark:     |
| DATE            | :white_check_mark:     |
| TIMESTAMP       | :white_check_mark:     |
| INTERVAL        | :white_check_mark:     |
| JSON            | :white_check_mark: (1) |
| XML             | :white_check_mark: (1) |
| MEM_TYPE        | :white_check_mark: (1) |
| CLOB            | :white_check_mark: (2) |
| NCLOB           | :white_check_mark: (2) |
| URI TYPE        | :white_check_mark: (2) |

* (1) supported as a string
* (2) you get unreadable representation, but in query you can cast such a types to varchar, hex

## Supported features

* Schema browsing
* Build query
* Visual query
* Connection test

## Remarks

* Do not finish your queries with ';' or you will get an exception specific to Oracle

See also:

* [Data connection](../data-connection.md)
* [Oracle Database](https://www.oracle.com/database/)
* [Oracle Database Wiki](https://en.wikipedia.org/wiki/Oracle_Database)
