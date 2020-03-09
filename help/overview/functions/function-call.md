<!-- TITLE: Function Call -->
<!-- SUBTITLE: -->

# Function Call

Function Call is a result of executing a [Data Job](../../access/data-job.md), [Data Query](../../access/data-query.md), 
[Script](../../compute/scripting.md) or any other [Function](function.md).

## Data

Each function call contains the following data:

  * Action
  * User that triggered job execution
  * Started on
  * Completed on
  * Status
  * [Runs](function-call.md) produced as a result of executing child actions

## Access control

Function Calls are first-class entities in the Datagrok platform, and as such are subjects to the 
standard checks and routines performed against them whenever they are used in the specific context. 
Some of the most popular privileges are:

  *  can_create
  *  can_edit
  *  can_delete
  *  can_query

Those privileges can be given to individuals or to groups (which can be defined via dynamic filters).
For more information on the access privilege model, refer to the Datagrok - Access Privileges page.

## Filtering

You can use these fields to filter action runs with [smart search](../smart-search.md):

| Field       | Description                                        |
|-------------|----------------------------------------------------|
| id          |                                                    |
| name        |                                                    |
| action      | [Func](function.md) object                      |
| childRuns   | list of [FuncCall](function-call.md) object      |
| parentRun   | [FuncCall](function-call.md) object              |
| status      |                                                    |
| started     |                                                    |
| finished    |                                                    | 
| createdOn   |                                                    |
| updatedOn   |                                                    | 


See also:

  * [Data Pipeline](../../access/data-pipeline.md)
  * [Data Source](../../access/data-source.md)
  * [Data Connection](../../access/data-connection.md)
  * [Data Query](../../access/data-query.md)
  * [Data Job](../../access/data-job.md)
