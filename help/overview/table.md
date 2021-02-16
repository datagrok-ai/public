<!-- TITLE: Table -->
<!-- SUBTITLE: -->

# Table

Table (also known as a data frame) is a two-dimensional structure with named columns
of different types. Supported types are: `string`, `bool`, `int`, `bigint`, `double`, `qnum` and `datetime`.
  
## Column

Dataframes consist of columns. In addition to the data type, a column can be associated with [tags](../discover/tags.md) that specify units, data format, semantic type, tooltip content and so on (right-click on a column header and open `Properties...`).
  
## Metadata

Each table can be annotated with attributes in the form of key-value pairs. Datagrok automatically annotates
certain tables with metadata (such as the source of the table, or time of import, etc). To edit metadata manually,
use 'Properties' context menu option.

It is possible to use metadata as a search criteria in the [Projects](../overview/project.md)

## Filtering

You can use these fields to filter tables with [smart search](smart-search.md):

| Field       | Description                                        |
|-------------|----------------------------------------------------|
| id          |                                                    |
| name        |                                                    |
| rowCount    |                                                    |
| colCount    |                                                    |
| createdOn   |                                                    |
| updatedOn   |                                                    | 
| author      | [User](../govern/user.md) object                             |
| starredBy   | [User](../govern/user.md) object                             |
| commentedBy | [User](../govern/user.md) object                             |


See also:

  * [JS API: DataFrame](https://datagrok.ai/js-api/DataFrame)
  * [Project](project.md)
  * [Table View](table-view.md)
  * [View Layout](../visualize/view-layout.md)
  * [Tags](../discover/tags.md)
  * [Grid](../visualize/viewers/grid.md)
