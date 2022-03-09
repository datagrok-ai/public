<!-- TITLE: Search patterns -->
<!-- SUBTITLE: -->

# Search patterns

Search patterns let you use a commonly accepted notation to specify conditions in free text. Use the same syntax to
query in-memory datasets via the [search mechanism](data-search.md), and to query external databases
with [parameterized queries](../access/parameterized-queries.md). When querying databases, behind the scenes the
platform will parse the free-text query, and then execute a parameterized, safe, provider-specific SQL query on the
backend.

For searching within the in-memory table, patterns can be preceded with column name
(example: "age > 21"). External query patterns do not allow the column name to be specified, as the column name is
already specified at the query level.

## Numerical patterns

| Pattern                | Example    |
|------------------------|------------|
| equals                 | 5          |
| equals                 | = 5        |
| not equals             | != 5       |
| greater than           | > 5        |
| greater than or equals | >= 5       |
| less than              | < 5        |
| less than or equals    | <= 5       |
| range (inclusive)      | 10-20      |
| in                     | in (5, 10) |
| not in                 | not in (5, 10) |

If the input does not match above-mentioned patterns, 'exact value' search is used.

## String patterns

| Pattern                | Example             |
|------------------------|---------------------|
| contains               | contains foo        |
| starts with            | starts with Charles |
| ends with              | ends with District  |
| regular expression     |  .*\[0-9\]+         |
| in                     | in (asian, other)   |
| not in                 | not in (m, f)       |

String matching is case-insensitive. If the input does not match above-mentioned patterns, 'exact value' search is used.

## Datetime patterns

| Pattern                | Example                                       |
|------------------------|-----------------------------------------------|
| 1984                   | years 1984 (ignore date and time)             |
| 1984-1986              | years 1984, 1985, 1986 (ignore date and time) |
| June 1984              | June 1984 (ignore date and time)              |
| Oct 17, 2019           | date = 10/17/2019 (ignore time)               |
| 10/17/2019 5:24 pm     | exact date and time                           |
| before 10/17/2019      | before the specified date                     |
| after 10/17/2019       | after the specified date                      |
| today                  |                                               |
| this week              |                                               |
| this month             |                                               |
| this year              |                                               |
| yesterday              |                                               |
| last week              |                                               |
| last month             |                                               |
| last year              |                                               |

## Provider compatibility

| Source                 | Support |
|------------------------|---------|
| Grok in-memory         | +       |
| Access                 | +       |
| Athena                 | +       |
| Cassandra              |         |
| DB2                    | +       |
| Firebird               | +       |
| HBase                  | +       |
| Hive                   | +       |
| Hive2                  | +       |
| Oracle                 | +       |
| MariaDB                | +       |
| MS SQL                 | +       |
| MongoDB                |         |
| MySql                  | +       |
| Postgres               | +       |
| SQLite                 | +       |
| Teradata               | +       |
| Vertica                | +       |
| Redis                  |         |
| SPARQL                 | +       |
| SAP                    |         |
| SAS                    |         |

See also:

* [Data search](data-search.md)
* [Parameterized queries](../access/parameterized-queries.md)
