<!-- TITLE: Sharing -->
<!-- SUBTITLE: -->

# Sharing

Many types of objects within the Datagrok platform can be shared with other users or [groups](../govern/group.md). Such
shareable objects are called [entities](../overview/objects.md). When an object is shared, you are essentially granting
a [privilege](../govern/authorization.md) (typically, 'view' or 'edit') to a grantee. See
the [Security](../govern/security.md) article for details on how to manage groups and privileges.

To share an object, either right-click on it and choose `Share...`, or expand the `Sharing` section in the object's
property panel and choose `Share...` there.

![Share Connection](share-connection.gif)

![Property Panel Share](property-panel-share.png)

## Sharing connections and queries

Access rights of a database connection inherit access rights of a query. However, access rights of the query don't
inherit access rights of the database connection. Thus, if one shares a query, the associated database connection shall
automatically be shared. At the same time, when you are sharing a database connection, your queries aren't going to be
shared automatically.

As for web queries, they are automatically shared along with sharing the corresponding connection.

## Videos

[![Sharing](../uploads/youtube/sharing.png "Open on Youtube")](https://www.youtube.com/watch?v=7MBXWzdC0-I&t=1334s)

See also:

* [Entities](../overview/objects.md)
* [Privileges](../govern/authorization.md)
* [Security](../govern/security.md)
