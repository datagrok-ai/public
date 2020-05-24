<!-- TITLE: Tags -->

# semantic-type

Applicable to columns. Determines column's [semantic type](semantic-types.md); used for
making automatic suggestions on functions, layouts, and predictive models
applicable to the current context. Also, see [units](#units)

# units

Applicable to columns. Determines units in which the values are stored. 
Also, see [semantic-type](#semantic-type).

# format

Applicable to numeric or datetime columns. Determines the way values should 
be represented. In addition to the standard notation (example: "#.000"), the 
following formats are predefined and could be used (they also appear as choices
when you right-click on a column and open "Format" section):
 
* int
* two digits after comma
* four digits after comma
* max two digits after comma
* scientific
* money
* compact long
* compact
* compact simple currency
* percent
* thousand separator
* full precision

# formula

Formula used for creating a derived column. Edit it in the "Formula" section on the property panel
to recalculate.

# layout_id

Used for matching layout columns with table columns when a layout is applied.
Applicable to columns.  
 
# .tooltip

Applicable to tables. Contains comma-separated list of column names to be used as a row tooltip.
This tag starts with ".", therefore it is not shown in the UI. 

# .row.group.tooltip

JSON-serialized settings of the viewer that is used to visualize a group of rows on a
tooltip. It is shown when user moves the cursor over the area that represents multiple
rows (such a a histogram bin, or a pie chart pie). 
Applicable to tables.

# .semantic-detection-duration

Indicates number of milliseconds spent for detecting column's semantic type.
This tag starts with ".", therefore it is not shown in the UI.
Applicable to columns.  

# .choices

JSON-encoded list of strings used by the grid cell editor to populate a combo box.
See also [auto-choices](#.auto-choices).
Applicable to string columns.

# .auto-choices

When set to 'true', switches the cell editor to a combo box that only allows to choose values
from a list of already existing values in the column.
See also [choices](#.choices).
Applicable for string columns.

# cell-renderer

Column cell renderer.

# query

A [query](../access/data-query.md) that was used to produce this table.
Applicable to tables.

# import-time

When was the table created.
Applicable to tables.

# db

Applicable to columns. A database the content of this column was retrieved from.
Used for data augmentation and impact analysis.
See also [db-schema](#db-schema), [db-table](#db-table), [db-column](#db-column), [db-path](#db-path). 

# db-schema

Applicable to columns. A database schema the content of this column was retrieved from.
Used for data augmentation and impact analysis.
See also [db](#db), [db-table](#db-table), [db-column](#db-column), [db-path](#db-path). 

# db-table

Applicable to columns. A database table the content of this column was retrieved from.
Used for data augmentation and impact analysis.
See also [db](#db), [db-schema](#db-schema), [db-column](#db-column), [db-path](#db-path). 

# db-column

Applicable to columns. A database column the content of this column was retrieved from.
Used for data augmentation and impact analysis.
See also [db](#db), [db-schema](#db-schema), [db-table](#db-table), [db-path](#db-path). 

# db-path

Applicable to columns. A database column, in the form on "db.schema.table.column", that 
the content of this column was retrieved from. Used for data augmentation and impact analysis.
See also [db](#db), [db-schema](#db-schema), [db-table](#db-table), [db-column](#db-column). 

# id

Entity id, as it is stored in the database.
Applicable to tables, columnms, and other entities.

# data-connection-id

Id of the [data connection](../access/data-connection.md) that was used to populate the table.
Applicable to tables. 


See also:
* [Metadata](metadata.md)