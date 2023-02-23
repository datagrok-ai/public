---
title: "Postgres"
---

Provides access to [PostgreSQL](https://www.postgresql.org/) database using SQL queries via JDBC driver .

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
| `list<string>`         |             |                            |

## Supported output type

| Type                 | Supported          |
|----------------------|--------------------|
| bigInt               | :white_check_mark: |
| integer              | :white_check_mark: |
| double precision     | :white_check_mark: |
| real                 | :white_check_mark: |
| numeric              | :white_check_mark: |
| serial               | :white_check_mark: |
| monetary             | not tested         |
| character type       | :white_check_mark: |
| bytea                | limited support    |
| date/time            | :white_check_mark: |
| boolean              | :white_check_mark: |
| network address type | :white_check_mark: |
| bit string           | :white_check_mark: |
| uuid                 | :white_check_mark: |
| jsonb, json          | limited support    |
| xml                  | limited support    |
| composite types      | limited support    |

## Supported features

* Schema browsing
* Connection test

See also:

* [Data connection](../data-connection.md)
* [PostgreSQL](https://www.postgresql.org/)
* [Wikipedia](https://en.wikipedia.org/wiki/PostgreSQL)
