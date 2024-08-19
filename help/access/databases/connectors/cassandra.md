---
title: "Cassandra"
---

Provides access to [Apache Cassandra](https://cassandra.apache.org/) database
using SQL queries via JDBC driver.

## Connection parameters

```json
{
  "server": "",
  "port": "",
  "db": "",
  "connString": ""
}
```

## Supported Parameters (1)

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
| `datetime`             | before      | before 1/1/2022            |
|                        | after       | after 1/1/2022             |
|                        | today       |                            |
|                        | this week   |                            |
|                        | this month  |                            |
|                        | this year   |                            |
|                        | last year   |                            |
|                        | min-max     |                            |
|                        | April 2021  |                            |

* (1) If you want to use parameterized queries, you need to add `ALLOW FILTERING`
  which can lead to sirius performance penalties. Moreover, to be able to use string patterns you need to create secondary indexes for columns in where clause.
  See [documentation](https://docs.datastax.com/en/cql-oss/3.3/cql/cql_reference/cqlSelect.html).

## Supported output types

| Type                                   | Supported              |
|----------------------------------------|------------------------|
| text                                   | :white_check_mark:     |
| tinyint, smallint, int, bigint, varint | :white_check_mark:     |
| decimal, float, double                 | :white_check_mark:     |
| boolean                                | :white_check_mark:     |
| date, timestamp, time                  | :white_check_mark:     |
| uuid                                   | :white_check_mark:     |
| duration                               | :white_check_mark: (1) |
| map, list, set, udt                    | :white_check_mark: (1) |
| blob                                   | limited support    (2) |

* (1) supported as a string
* (2) you get unreadable representation, but in query you can cast such a types to varchar, hex

## Supported features

* Schema browsing
* Connection test

See also:

* [Data connection](../../access.md#data-connection)
* [Apache Cassandra](https://cassandra.apache.org/)
