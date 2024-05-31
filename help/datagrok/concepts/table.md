---
title: "Dataframe"
sidebar_position: 1
---

Table (also known as a dataframe) is a two-dimensional structure with named columns of different types. Supported types
are: `string`, `bool`, `int`, `bigint`, `double`, `qnum`, `datetime`.

## Column

Dataframes consist of columns. In addition to the data type, a column can be
associated with [tags](../../govern/catalog/tags.md) that specify units, data
format, [semantic type](../../govern/catalog/semantic-types.md), tooltip content
and so forth. 

To view column properties, right-click a column header and select
**Properties...** The [Context Panel](../navigation/panels/panels.md#context-panel) updates to shows the column's properties.

## Metadata

Like columns, each table can be annotated with attributes in the form of
key-value pairs. Datagrok automatically annotates certain tables with metadata
(such as the source of the table, or time of import, etc). To edit metadata
manually, right-click a column header and select
**Column Properties...**.

<!--TODO: revise when this feature (manually add metadata) is updated. Possibly merge Column Properties and Properties context menu options -->

You can search tables by metadata in [projects](project/project.md).

## Filtering

You can use these fields to filter tables with smart search:

| Field       | Description                      |
|-------------|----------------------------------|
| ID          |                                  |
| name        |                                  |
| rowCount    |                                  |
| colCount    |                                  |
| createdOn   |                                  |
| updatedOn   |                                  |
| author      | [User](../../govern/access-control/users-and-groups#users) object |
| starredBy   | [User](../../govern/access-control/users-and-groups#users) object |
| commentedBy | [User](../../govern/access-control/users-and-groups#users) object |

See also:

* [JS API: Dataframe](https://datagrok.ai/js-api/classes/dg.DataFrame)
* [Project](project/project.md)
* [Table view](../navigation/views/table-view.md)
* [View layout](../../visualize/view-layout.md)
* [Tags](../../govern/catalog/tags.md)
* [Grid](../../visualize/viewers/grid.md)
