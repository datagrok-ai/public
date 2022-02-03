<!-- TITLE: Tags -->

# Tags

Most of the objects in Datagrok can be annotated with metadata (key-value pairs). The metadata could be set manually;
additionally, some of it gets assigned automatically. Some keys affect the way an object (such as a column) interacts
with the platform; other have no effect at all, except that you can search objects by metadata.

Below are some of the standard tags related to tables or columns. To edit column's metadata, right-click on it and
select "Properties..." (or press F2 in the grid).

## Quality

Applicable to columns. Determines column's [semantic type](semantic-types.md); used for making automatic suggestions on
functions, layouts, and predictive models applicable to the current context. Also, see [units](#units).

## Units

Applicable to columns. Determines units in which the values are stored. Also, see [quality](#quality).

## Format

Applicable to numeric or datetime columns. Determines the way values should be represented.

### Numbers

In addition to the standard notation (example: "#.000"), the following formats are predefined and could be used (they
also appear as choices when you right-click on a column and open "Format"
section):

|           Format           | Example |
|----------------------------|---------|
| int                        | 71      |
| two digits after comma     | 70.50   |
| four digits after comma    | 70.5000 |
| max two digits after comma | 70.5    |
| scientific                 | 7E1     |
| money                      | $70.50  |
| compact long               | 70.5    |
| compact                    | 70.5    |
| compact simple currency    | $70.50  |
| percent                    | 7,050%  |
| thousand separator         | 71      |
| full precision             | 70.5    |

### Datetime

The standard formats representing date and time include:

|         Format          |          Example            |
|-------------------------|-----------------------------|
| MM/dd/yyyy HH:mm:ss.fff | 07/21/1991 00:00:00.000     |
| M/d/yyyy                | 7/21/1991                   |
| dd.MM.yyyy              | 21.07.1991                  |
| M/d/yyyy h:mm tt        | 7/21/1991 12:00 AM          |
| M/d/yyyy h:mm:ss tt     | 7/21/1991 12:00:00 AM       |
| yyyy-MM-dd              | 1991-07-21                  |
| dddd, MMMM d, yyyy      | Sunday, July 21, 1991       |
| MMM d, yyyy             | Jul 21, 1991                |
| h:mm:ss                 | 12:00:00                    |
| h:mm:ss.fff             | 12:00:00.000                |
| relative                | 29 years ago                |
| auto                    | Jul 21, 1991                |

In case you need to define a different date format, proceed to column properties and edit the template string. Here are
the valid specifiers:

| Symbol |          Meaning                 |          Example           |
|--------|----------------------------------|----------------------------|
| yy     | Year without the century         | 00, 01, ..., 20, ..., 99   |
| yyyy   | Year with the century            | 0001, ..., 2020, ..., 9999 |
| M      | Month                            | 1, 2, 3, ..., 12           |
| MM     | Zero-padded month                | 01, 02, 03, ..., 12        |
| MMM    | Abbreviated month name           | Jan, Feb, Mar, ..., Dec    |
| MMMM   | Full month name                  | January, ..., December     |
| d      | Day of the month                 | 1, 2, 3, ..., 31           |
| dd     | Zero-padded day                  | 01, 02, 03, ..., 31        |
| ddd    | Abbreviated weekday name         | Mon, ..., Fri, Sat, Sun    |
| dddd   | Full weekday name                | Monday, ..., Sunday        |
| h      | Hour (12-hour clock)             | 1, 2, 3, ..., 12           |
| hh     | Zero-padded hour (12-hour clock) | 01, 02, 03, ..., 12        |
| H      | Hour (24-hour clock)             | 0, 1, 2, ..., 23           |
| HH     | Zero-padded hour (24-hour clock) | 00, 01, 02, ..., 23        |
| m      | Minute                           | 0, 1, 2, ..., 59           |
| mm     | Zero-padded minute               | 00, 01, 02, ..., 59        |
| s      | Second                           | 0, 1, 2, ..., 59           |
| ss     | Zero-padded second               | 00, 01, 02, ..., 59        |
| f      | Second fraction (1 digit)        | 0, 1, ..., 9               |
| ff     | Second fraction (2 digits)       | 00, 01, ..., 99            |
| fff    | Second fraction (3 digits)       | 000, 001, ..., 999         |
| tt     | 12-hour periods                  | AM, PM                     |

## .color-coding-type

Indicates on what basis the columns are colored. Numerical columns have the options `Off`,
`Linear`, and `Conditional`, while categorical columns are limited to `Off` and `Categorical`.

## .color-coding-conditional

Applies conditional formatting to numeric and categorical columns. See
an [example](https://public.datagrok.ai/js/samples/grid/color-coding-conditional).

## Formula

Formula used for creating a derived column. Edit it in the "Formula" section on the property panel to recalculate. Note
that changing the `formula` tag does not cause recalculation.

## Layout-id

Applicable to columns. Used for matching layout columns with table columns when a layout is applied.
See [layout suggestions](../visualize/view-layout.md#layout-suggestions) for details.

## .tooltip

Applicable to tables. Contains comma-separated list of column names to be used as a row tooltip. This tag starts with "
.", therefore it is not shown in the UI.

## .row.group.tooltip

JSON-serialized settings of the viewer that is used to visualize a group of rows on a tooltip. It is shown when user
moves the cursor over the area that represents multiple rows (such a a histogram bin, or a pie chart pie). Applicable to
tables.

## .semantic-detection-duration

Indicates number of milliseconds spent for detecting column's semantic type. This tag starts with "
.", therefore it is not shown in the UI. Applicable to columns.

## .choices

JSON-encoded list of strings used by the grid cell editor to populate a combo box. See
also [auto-choices](#.auto-choices). Applicable to string columns.

## .auto-choices

When set to 'true', switches the cell editor to a combo box that only allows to choose values from a list of already
existing values in the column. See also [choices](#.choices). Applicable for string columns.

## Cell-renderer

Column cell renderer.

## Query

A [query](../access/data-query.md) that was used to produce this table. Applicable to tables.

## Import-time

When was the table created. Applicable to tables.

## Db

Applicable to columns. A database the content of this column was retrieved from. Used for data augmentation and impact
analysis. See also [db-schema](#db-schema), [db-table](#db-table)
, [db-column](#db-column), [db-path](#db-path).

## Db-schema

Applicable to columns. A database schema the content of this column was retrieved from. Used for data augmentation and
impact analysis. See also [db](#db), [db-table](#db-table)
, [db-column](#db-column), [db-path](#db-path).

## Db-table

Applicable to columns. A database table the content of this column was retrieved from. Used for data augmentation and
impact analysis. See also [db](#db), [db-schema](#db-schema)
, [db-column](#db-column), [db-path](#db-path).

## Db-column

Applicable to columns. A database column the content of this column was retrieved from. Used for data augmentation and
impact analysis. See also [db](#db), [db-schema](#db-schema)
, [db-table](#db-table), [db-path](#db-path).

## Db-path

Applicable to columns. A database column, in the form on "db.schema.table.column", that the content of this column was
retrieved from. Used for data augmentation and impact analysis. See also [db](#db), [db-schema](#db-schema)
, [db-table](#db-table), [db-column](#db-column).

## ID

Entity id, as it is stored in the database. Applicable to tables, columnms, and other entities.

## Data-connection-id

ID of the [data connection](../access/data-connection.md) that was used to populate the table. Applicable to tables.

## .history

History of all modifications applied to that table. Applies to tables.

## .script

[Grok script](../overview/grok-script.md) that was used to create a table. It could represent getting a data via
a [database query](../access/data-query.md), from a web service, from a file share, or using any
other [function](../overview/functions/function.md) that returns a [table](../overview/table.md).

If this tag is present in a table, a "data sync" option appears next to this table in the
"Upload table" dialog. If the option is checked, the table will not be uploaded to the server. Instead, the script will
get re-executed when user opens this project next time.

Applies to tables.

## Chem-descriptor

A [molecular descriptor](../domains/chem/descriptors.md) used for calculating the values of that column.

## Chem-fingerprinter

A [molecular fingerprinter](../domains/chem/fingerprints.md) used for calculating the values of that column.

See also:

* [Metadata](metadata.md)
* [JS API: metadata](https://public.datagrok.ai/js/samples/data-frame/metadata)
