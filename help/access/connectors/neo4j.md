---
title: "Neo4j"
---

This is a [connector](../data-connection.md#connectors) that provides access to the [Neo4j](https://neo4j.com/) graph
database via JDBC driver. Allows to query Neo4j using [Cypher](https://neo4j.com/developer/cypher-query-language)
language, and use results in dashboards, data augmentation panels, or via the [JS API](../../develop/js-api.md).

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

* Note, that patterns are supported only in WHERE clause, not in MATCH

## Supported output property types

| Type                    | Supported              |
|-------------------------|------------------------|
| Boolean                 | :white_check_mark:     |
| Date                    | :white_check_mark:     |
| DateTime, LocalDateTime | :white_check_mark:     |
| LocalTime, Time         | :white_check_mark:     |
| Float                   | :white_check_mark:     |
| Integer                 | :white_check_mark:     |
| String                  | :white_check_mark:     |
| Point                   | :white_check_mark: (1) |

* (1) supported as a string

## Supported features

* Connection test

## Remarks

* There is one limitation due to JDBC driver of Neo4j. If your query return 0 rows you will get an exception.
* Neo4j uses ``//`` for comments. Use this for parameters declaration in query too.

See also:

* [Data connection](../data-connection.md)
* [Neo4j](https://neo4j.com/)
