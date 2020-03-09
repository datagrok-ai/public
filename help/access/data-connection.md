<!-- TITLE: Data Connection -->
<!-- SUBTITLE: -->

# Data Connection

Data connection defines a way to connect to a particular data source. Connection parameters depend on 
the data source. Typically, you would need to provide server name and login credentials.

Connections are first-class entities in the Datagrok platform, and as such are subjects to the standard 
checks and routines performed against them whenever they are used in the specific context. Some of the 
most popular privileges are:

  * can_view
  * can_edit
  * can_delete
  * can_share


Those privileges can be given to individuals or to groups (which can be defined via dynamic filters). For more 
information on the access privilege model, refer to the Datagrok - Access Privileges page.

Another “out of the box” feature that comes with connections being first-class entity is the audit trail for 
every action performed against the connection. For details on that, check out [Audit](../govern/audit.md) page.

## Filtering

You can use these fields to filter connections with [smart search](../overview/smart-search.md):

| Field       | Description                                 |
|-------------|---------------------------------------------|
| id          |                                             |
| name        |                                             |
| server      |                                             |
| port        |                                             |
| db          |                                             |
| login       |                                             |
| dataSource  |                                             |
| description |                                             |
| createdOn   |                                             |
| updatedOn   |                                             |
| author      | [User](../govern/user.md) object                      |
| starredBy   | [User](../govern/user.md) object                      |
| commentedBy | [User](../govern/user.md) object                      |
| usedBy      | [User](../govern/user.md) object                      |


See also:

  * [Data Pipeline](data-pipeline.md)
  * [Data Source](data-source.md)
  * [Data Query](data-query.md)
  * [Data Job](data-job.md)
  * [Function Call](../overview/functions/function-call.md)
