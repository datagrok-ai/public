---
title: "Vertica"
---

Provides access to [Vertica](https://www.vertica.com/overview/) database using
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
|                        | \<=          | \<=100                      |
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

| Type                                                                                   | Supported              |
|----------------------------------------------------------------------------------------|------------------------|
| CHAR, VARCHAR, LONG VARCHAR                                                            | :white_check_mark:     |
| APPROXIMATE NUMERIC, EXACT NUMERIC                                                     | :white_check_mark:     |
| DATE, TIME, TIMESTAMP, TIME <br/>WITH TIMEZONE, <br/>TIMESTAMP WITH TIMEZONE, INTERVAL | :white_check_mark:     |
| BOOLEAN                                                                                | :white_check_mark:     |
| UUID                                                                                   | :white_check_mark:     |
| SPATIAL                                                                                | limited support    (1) |
| BINARY, VARBINARY, LONG VARBINARY                                                      | limited support    (1) |

* (1) you get unreadable representation, but in query you can cast such a types to varchar, hex

## Supported features

* Schema browsing
* Join DB tables
* Aggregation query
* Connection test

See also:

* [Data connection](../../access.md#data-connection)
* [Vertica](https://www.vertica.com/overview/)
