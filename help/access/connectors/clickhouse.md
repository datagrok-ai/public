---
title: "ClickHouse"
---

Provides access to the [ClickHouse](https://clickhouse.com/clickhouse) database using SQL queries via JDBC driver.

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

| Type                                            | Supported              |
|-------------------------------------------------|------------------------|
| Int8, Int16, Int32, Int64, Int128, Int256       | :white_check_mark:     |
| UInt8, UInt16, UInt32, UInt64, UInt128, UInt256 | :white_check_mark:     |
| Float32, Float64                                | :white_check_mark:     |
| Decimal32, Decimal64, Decimal128, Decimal256    | :white_check_mark:     |
| Bool                                            | :white_check_mark:     |
| String                                          | :white_check_mark:     |
| UUID                                            | :white_check_mark:     |
| Date, Date32, DateTime, DateTime64              | :white_check_mark:     |
| Tuple, Map, Array, Nested           (1)         | :white_check_mark: (1) |
| Point, Ring, Polygon, MultiPolygon  (1)         | :white_check_mark: (1) |

* (1) supported as a string

## Supported features

* Schema browsing
* Build query
* Visual query
* Connection test

See also:

* [Data connection](../data-connection.md)
* [ClickHouse](https://clickhouse.com/clickhouse)
