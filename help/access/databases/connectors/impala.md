---
title: "Impala"
---

Provides access to [Impala](https://impala.apache.org/) query engine via JDBC driver.

## Connection parameters

```json
{
  "server": "",
  "port": "",
  "schema": "",
  "connString": ""
}
```

## Supported Parameters

| Type                   | Value                | Description or Example     |
|------------------------|----------------------|----------------------------|
| `num`, `int`, `double` | =                    | =100                       |
|                        | >                    | >1.02                      |
|                        | >=                   | >=4.1                      |
|                        | \<=                   | \<=100                      |
|                        | !=                   | !=5                        |
|                        | in                   | in (1, 3, 10.2)            |
|                        | min-max              | 1.5-10.0                   |
|                        | is null/ is not null |                            |
| `string`               | contains             | contains ea                |
|                        | starts with          | starts with R              |
|                        | ends with            | ends with w                |
|                        | in                   | in (ab, "c d", "e\\"f\\"") |
|                        | regex                | regex ^(.+)@(.+)$          |
|                        | is null/ is not null |                            |
| `datetime`             | anytime              |                            |
|                        | before               | before 1/1/2022            |
|                        | after                | after 1/1/2022             |
|                        | today                |                            |
|                        | this week            |                            |
|                        | this month           |                            |
|                        | this year            |                            |
|                        | last year            |                            |
|                        | min-max              |                            |
|                        | April 2021           |                            |
|                        | last                 | last 10 days, last 2 weeks |
|                        | is null/ is not null |                            |
| `list<string>` (1)     |                      |                            |

* (1) default parameters are not supported

## Supported output types

| Type                           | Supported              |
|--------------------------------|------------------------|
| BIGINT, INT, SMALLINT, TINYINT | :white_check_mark:     |
| BOOLEAN                        | :white_check_mark:     |
| CHAR, VARCHAR, STRING          | :white_check_mark:     |
| DECIMAL, DOUBLE, FLOAT         | :white_check_mark:     |
| DATE, TIMESTAMP                | :white_check_mark:     |
| BINARY                         | limited support    (1) |

* (1) you get unreadable representation, but in query you can cast such a types to varchar, hex

## Supported features

* Schema browsing
* Join DB tables
* Aggregation query
* Connection test

See also:

* [Data connection](../../access.md#data-connection)
* [Impala](https://impala.apache.org/)
