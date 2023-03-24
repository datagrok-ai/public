---
title: "Snowflake"
---

Provides access to [Snowflake](https://www.snowflake.com/) database using SQL
queries via a JDBC driver.

## Connection parameters

```json
{
  "region": "",
  "accountLocator": "",
  "cloud": "",
  "db": "",
  "warehouse": "",
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
|                        | min-max     | 2021-2023                  |
|                        | April 2021  |                            |
| `list<string>` (1)     |             |                            |

* (1) default parameters are not supported

## Supported output types

| Type                | Supported              |
|---------------------|------------------------|
| NUMBER              | :white_check_mark:     |
| FLOAT               | :white_check_mark:     |
| VARCHAR             | :white_check_mark:     |
| BOOLEAN             | :white_check_mark:     |
| DATE                | :white_check_mark:     |
| TIMESTAMP           | :white_check_mark:     |
| VARIANT             | :white_check_mark: (1) |
| OBJECT              | :white_check_mark: (1) |
| ARRAY               | :white_check_mark: (1) |
| GEOGRAPHY           | :white_check_mark: (1) |
| GEOMETRY            | :white_check_mark: (1) |
| BINARY, VARBINARY   | limited support    (2) |

* (1) supported as a string
* (2) you get unreadable representation, but in query you can cast such a types to varchar, hex

## Supported features

* Schema browsing
* Build query
* Visual query
* Connection test

See also:

* [Data connection](../data-connection.md)
* [Snowflake](https://www.snowflake.com/)
