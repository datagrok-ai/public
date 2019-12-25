<!-- TITLE: Table -->
<!-- SUBTITLE: -->

# Table

Table (also known as a data frame) is a two-dimensional structure with named columns
of different types. Supported types are: string, int, float, bool, DateTime.
  
## Column

TODO: GROK-3335
  
## Metadata

Each table can be annotated with attributes in the form of key-value pairs. Datagrok automatically annotates
certain tables with metadata (such as the source of the table, or time of import, etc). To edit metadata manually,
use 'Properties' context menu option.

It is possible to use metadata as a search criteria in the [Projects](../views/welcome-view.md#Projects)

## Filtering

You can use these fields to filter tables with [smart search](../features/smart-search.md):

| Field       | Description                                        |
|-------------|----------------------------------------------------|
| id          |                                                    |
| name        |                                                    |
| rowCount    |                                                    |
| colCount    |                                                    |
| createdOn   |                                                    |
| updatedOn   |                                                    | 
| author      | [User](user.md) object                             |
| starredBy   | [User](user.md) object                             |
| commentedBy | [User](user.md) object                             |


See also:

  * [Project](project.md)
  * [Table View](../views/table-view.md)
  * [View Layout](../entities/view-layout.md)
