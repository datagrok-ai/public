<!-- TITLE: Tests: Build Query -->
<!-- SUBTITLE: -->

# Tests: Data query

[Query Builder](../../access/query-builder.md) allows to build query using visual interface.

## Testing scenario

1. Open "Connect to data..." tab on "Welcome" view

1. Expand ```PosgteSQL -> northwind -> Tables```

1. Open [Build Query](../../access/query-builder.md) dialog for **employees** DB table (from it's context menu)

* [Build Query](../../access/query-builder.md) dialog is open
* Context help switched to [Build Query](../../access/query-builder.md)

1. Select all columns from **employees** table for query

* SQL query added to dialog field (*)
* Result preview is shown in dialog

1. Click on "Save as query" from dialog menu (v)

* Added "Query View", in which query from [Query Builder](../../access/query-builder.md)
* "Name" field matches DB table name (*employees*)

1. Click on "Add result to workspace" from dialog menu (v)

* Result table of query from [Query Builder](../../access/query-builder.md) has been added to workspace

1. Return to view with [Source Tree](../../access/data-source.md)

1. Click on "Get All" from context menu of **employees** DB table

* Table "*employees*" with all values and columns added to workspace

1. Repeat previous steps for **Postgres**, **MySql**, **MS SQL**, **MariaDB**, **ORACLE**
   providers

(*):

```
select
  employees.address,
  employees.birthdate,
  employees.city,
  employees.country,
  employees.employeeid,
  employees.extension,
  employees.firstname,
  employees.hiredate,
  employees.homephone,
  employees.lastname,
  employees.notes,
  employees.photo,
  employees.photopath,
  employees.postalcode,
  employees.region,
  employees.reportsto,
  employees.title,
  employees.titleofcourtesy
from
  employees
```

See also:

* [Data Sourse test](../../access/data-source-test.md)
* [Data Query](../../access/data-query.md)
