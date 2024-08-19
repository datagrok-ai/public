<!-- TITLE: Tests: Join DB Tables -->
<!-- SUBTITLE: -->

# Tests: Data query

[Query Builder](../../access/databases/databases.md#join-tables) allows to build query using visual interface.

## Testing scenario

1. Open "Connect to data..." tab on "Welcome" view

1. Expand ```PosgteSQL -> northwind -> Tables```

1. Open [Join DB Tables](../../access/databases/databases.md#join-tables) dialog for **employees** DB table (from it's context menu)

* [Join DB Tables](../../access/databases/databases.md#join-tables) dialog is open
* Context help switched to [Join DB Tables](../../access/databases/databases.md#join-tables)

1. Select all columns from **employees** table for query

* SQL query added to dialog field (*)
* Result preview is shown in dialog

1. Click on "Save as query" from dialog menu (v)

* Added "Query View", in which query from [Query Builder](../../access/databases/databases.md#join-tables)
* "Name" field matches DB table name (*employees*)

1. Click on "Add result to workspace" from dialog menu (v)

* Result table of query from [Query Builder](../../access/databases/databases.md#join-tables) has been added to workspace

1. Return to view with [Source Tree](../../access/access.md#data-sources)

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
* [Data Query](../../access/access.md#data-query)
