---
title: "Entities"
sidebar_position: 3
---

Certain classes of objects in Datagrok have a common set of operations and
features applicable to them. We call these objects _entities_. Here they are:

* [Function](functions/functions.md)
* [Function Call](functions/function-call.md)
* Data Pipeline
* [Data Connection](../../access/access.md#data-connection)
* [Data Query](../../access/access.md#data-query)
* Data Job
* DB Table Info
* DB Column Info
* [User](../../govern/access-control/users-and-groups.md#users)
* [Group](../../govern/access-control/users-and-groups.md#groups)
* [Model](../../learn/learn.md)
* [Notebook](../../compute/jupyter-notebook.md)
* [Package](../../develop/develop.md#packages)
* [Project](project/project.md)
* Repository
* [Script](../../compute/scripting/scripting.mdx)
* [Table (dataframe)](table.md)
* [View Layout](../../visualize/view-layout.md)

The following operations can be applied to any entity:

* Getting its URL
* Referencing it in a [chat](../../collaborate/chat.md), in
  [markup](../../develop/under-the-hood/markup.md), or a [dashboard](project/dashboard.md)
* Assigning [privileges](../../govern/access-control/access-control.md#authorization), such as rights to view,
  edit, or share, to a particular instance
* Using it as a parameter in the [audit](../../govern/audit/audit.md) record
* Deleting it

## Metadata

Any entity can be annotated with metadata, which you can use to search entities. In
Datagrok, there are three kinds of metadata:

* Tags
* Parameters
* Properties

### Tags

A tag is a short string that can be associated with any entity. Use it to
organize and group different entities. For example, a `#chem` tag could
be used with chemical scripts, molecular datasets, or chemical database queries.

To search for tagged entities, enter `#tagname` in the search box.

### Parameters

Each entity can be associated with a list of key-value pairs.

To search for entities with the specified parameters, use the `<paramName> <op> <value>` syntax, where `op` is an operator like `>` `!=`.

For example:

* `imported < 10/02/2019`
* `rows > 200`
* `author = @current`

### Properties

Properties are similar to parameters, but unlike parameters, properties can
include nested properties. Additionally, you can't edit them directly. 

Some properties are universal across all entities. Others are specific to
particular entity types. Common properties include:

* ID
* name
* author
* createdOn
* updatedOn
* commentedBy
* usedBy
* starredBy

## Search

You can use metadata to [search for entities](../navigation/views/table-view.md#search)

See also:

* [Privileges](../../govern/access-control/access-control.md#authentication)
* [Sharing](../navigation/basic-tasks/basic-tasks.md#share)
* [Audit](../../govern/audit/audit.md)
