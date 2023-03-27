---
title: "Hive2"
---

Provides access to [Apache Hive](https://hive.apache.org/) NoSQL database using
SQL queries via HiveSever2 JDBC driver.

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

| Type                                 | Supported              |
|--------------------------------------|------------------------|
| TINYINT, SMALLINT                    | :white_check_mark:     |
| INT, BIGINT                          | :white_check_mark:     |
| FLOAT, DOUBLE, DECIMAL               | :white_check_mark:     |
| TIMESTAMP, DATE, INTERVAL            | :white_check_mark:     |
| STRING, VARCHAR, CHAR                | :white_check_mark:     |
| BOOLEAN                              | :white_check_mark:     |
| ARRAY<data_type>                     | :white_check_mark: (1) |
| MAP<primitive_type, data_type>       | :white_check_mark: (1) |
| STRUCT<col_name : data_type  ...>    | :white_check_mark: (1) |
| UNIONTYPE<data_type, data_type, ...> | not tested             |
| BINARY                               | not tested             |

* (1) supported as a string

## Supported features

* Schema browsing (you need credentials from metastore database)
* Connection test

See also:

* [Data connection](../data-connection.md)
* [Apache Hive](https://hive.apache.org/)
* [Apache Hive Wiki](https://en.wikipedia.org/wiki/Apache_Hive)
