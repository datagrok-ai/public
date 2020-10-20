<!-- TITLE: Sharing -->
<!-- SUBTITLE: -->

# Sharing

Many types of objects within the Datagrok platform can be shared with other users or groups. When an object is
shared, you are essentially granting a privilege (typically, 'view' or 'edit') to grantee. See
 [security](../govern/security.md) for details on how to manage groups and privileges.

## List of shareable objects

* [Project](../overview/project.md)
* [Data Connection](../access/data-connection.md)
* [Data Query](../access/data-query.md)
* [Data Job](../access/data-job.md)
* [Function Call](../overview/functions/function-call.md)

## Sharing connections and queries

Access rights of the database connection inherit access rights of the query. However, access rights of the query don't inherit access rights of the database connection. Thus, if one shares a query, the associated database connection shall automatically be shared. This has to do with the fact However, by sharing the database connection, the queries aren't going to be shared automatically.

As for web queries, they are automatically shared along with sharing the corresponding connection.

See also:

* [Privileges](../govern/authorization.md)