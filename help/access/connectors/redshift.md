---
title: "Redshift"
---

Provides access to [Amazon Redshift](https://aws.amazon.com/en/redshift/)
database using SQL queries via JDBC driver.

## Connection parameters

```json
{
  "server": "",
  "port": "",
  "db": "",
  "connString": ""
}
```

* Server parameter should be in a form of: ``examplecluster.abc123xyz789.us-west-2.redshift.amazonaws.com``
* Also note that credentials should be from admin user account of Redshift service

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

| Type                  | Supported              |
|-----------------------|------------------------|
| bigInt                | :white_check_mark:     |
| integer               | :white_check_mark:     |
| smallint              | :white_check_mark:     |
| double precision      | :white_check_mark:     |
| real                  | :white_check_mark:     |
| decimal               | :white_check_mark:     |
| char, varchar         | :white_check_mark:     |
| date, timestamp, time | :white_check_mark:     |
| boolean               | :white_check_mark:     |
| super                 | :white_check_mark: (1) |
| geometry / geography  | not tested             |
| varbyte               | not tested             |

* (1) supported as a string

## Supported features

* Schema browsing
* Build query
* Visual query
* Connection test

See also:

* [Data connection](../data-connection.md)
* [Amazon Redshift](https://aws.amazon.com/en/redshift/)
