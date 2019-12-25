<!-- TITLE: Data Query -->
<!-- SUBTITLE: -->

# Data Query

Data query defines which data should be extracted from the data source. For databases, that would 
typically be the SQL query; for Excel file, that would be sheet name, etc.

Queries can be executed either manually, or as part of the [Data Job](data-job.md). The result of 
executing a query is represented by the [Function Call](function-call.md).

## Parameterized queries

Queries can be parameterized. In the following example, the we are introducing 'employeeId' parameter: 

```$sql
--input: int employeeId = 5
SELECT * FROM Orders WHERE (employeeId = @employeeId)
```

For more details, see [Parameterized queries](connect/parameterized-queries.md) section.

## Access control

Data queries are first-class entities in the Datagrok platform, and as such are subjects to the 
standard checks and routines performed against them whenever they are used in the specific context. 
Some of the most popular privileges are:

  * can_create
  * can_edit
  * can_delete
  * can_query

Those privileges can be given to individuals or to groups (which can be defined via dynamic filters). 
For more information on the access privilege model, refer to the Datagrok - Access Privileges page.

## Filtering

You can use these fields to filter queries with [smart search](../features/smart-search.md):

| Field       | Description                                        |
|-------------|----------------------------------------------------|
| id          |                                                    |
| name        |                                                    |
| query       |                                                    |
| runs        | list of [FuncCall](function-call.md) object      |
| connection  | [DataConnection](data-connection.md) object        |
| jobs        | [DataJob](data-job.md) object                      |
| createdOn   |                                                    |
| updatedOn   |                                                    | 
| author      | [User](user.md) object                             |
| starredBy   | [User](user.md) object                             |
| commentedBy | [User](user.md) object                             |
| usedBy      | [User](user.md) object                             |

See also:

  * [Edit Data Query](../dialogs/edit-query.md)
  * [Data Pipeline](data-pipeline.md)
  * [Data Source](data-source.md)
  * [Data Connection](data-connection.md)
  * [Data Job](data-job.md)
  * [Function Call](function-call.md)
  * [Parameterized queries](connect/parameterized-queries.md)
  * [Query Builder](../dialogs/query-builder.md)
