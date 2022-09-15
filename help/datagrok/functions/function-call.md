<!-- TITLE: Function call -->
<!-- SUBTITLE: -->

# Function call

Function Call is a result of executing a [Data Job](../../access/data-job.md)
, [Data Query](../../access/data-query.md),
[Script](../../compute/scripting.md) or any other [Function](function.md).

## Data

Each function call contains the following data:

* Function
* User that triggered job execution
* Started on
* Completed on
* Status
* [Runs](function-call.md) produced as a result of executing child actions

## Access control

Connections are first-class entities in the Datagrok platform, and as such are subjects to the standard checks and
routines performed against them whenever they are used in the specific context. Some of the most popular privileges
are: `view`, `edit`, `delete`, and `share`. Those privileges can be given to individual users, or
to [groups](../../govern/group.md). For more information on the access privilege model, check
out [privileges](../../govern/security.md#privileges).

Those privileges can be given to individuals or to groups (which can be defined via dynamic filters)
. For more information on the access privilege model, refer to the Datagrok - Access Privileges page.

## Filtering

You can use these fields to filter action runs with [smart search](../smart-search.md):

| Field       | Description                                        |
|-------------|----------------------------------------------------|
| ID          |                                                    |
| name        |                                                    |
| action      | [Func](function.md) object                         |
| childRuns   | list of [FuncCall](function-call.md) object        |
| parentRun   | [FuncCall](function-call.md) object                |
| status      |                                                    |
| started     |                                                    |
| finished    |                                                    |
| createdOn   |                                                    |
| updatedOn   |                                                    |

See also:

* [Data pipeline](../../access/data-pipeline.md)
* [Data connection](../../access/data-connection.md)
* [Data query](../../access/data-query.md)
* [Data job](../../access/data-job.md)
