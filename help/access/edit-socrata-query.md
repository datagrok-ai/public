<!-- TITLE: Edit Socrata Query -->
<!-- SUBTITLE: -->

# Edit Socrata Query

Edit attributes of data SoQL query.

SoQL statements are broken into "parameters" similar to clauses in SQL statements. 
Each clause can be expressed either directly as a URL parameter or as a SoQL statement. 
If a parameter is not specified, then the default is used. 
See [Queries using SODA](https://dev.socrata.com/docs/queries) for more details:

| Parameter | Description | Default | In query |
|-----------|-------------|---------|-----------|
| select    | The set of columns to be returned, similar to a SELECT in SQL | All columns, equivalent to select=\* | SELECT |
| where     | Filters the rows to be returned, similar to WHERE | No filter | WHERE |
| order     | Column to order results on, similar to ORDER BY in SQL | Unspecified order | ORDER BY |
| group     | Column to group results on, similar to GROUP BY in SQL | No grouping | GROUP BY |
| having    | Filters the rows that result from an aggregation, similar to HAVING | No filter | HAVING |
| limit     | Maximum number of results to return | 1000 (2.0 endpoints: maximum of 50,000; 2.1: unlimited ) | LIMIT |
| offset    | Offset count into the results to start at, used for paging | 0 | OFFSET |
| q         | Performs a full text search for a value. | No search | N/A |
| query     | A full SoQL query string, all as one parameter | N/A | N/A |
| bom       | Prepends a UTF-8 Byte Order Mark to the beginning of CSV output | false | N/A |

Note that for equality comparisons, the $where clause can be replaced with using the column 
name as the query parameter.

See also:

  * [Data Query](data-query.md)
  * [Queries using SODA](https://dev.socrata.com/docs/queries)
  * [Data Pipeline](data-pipeline.md)
  * [Data Source](data-source.md)
  * [Data Connection](data-connection.md)
  * [Data Query](data-query.md)
  * [Data Job](data-job.md)
  * [Function Call](../overview/functions/function-call.md)
