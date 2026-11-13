# db-explorer changelog

## v.next

* Added catalogs: a configuration can name the catalog of its schema (`catalogName`) for a database that has them (Databricks, Snowflake, SQL Server); its tables, joins and references are then read as `catalog.schema.table`, and the editor offers the catalog before the schema and reads the schemas and tables of that catalog alone. A configuration without a catalog reads the connection's own catalog, as before
* Added binding the editor to one connection (`new DBExplorerEditor(onSave, false, connection)`): it then lists no other connection of the server, says which catalog and schema its tables are read from, and Import JSON brings a configuration made elsewhere, such as for the same tables in another catalog, into the schema it is open on
* Fixed the editor's preview cards failing on SQL Server (`limit 1`), and a quote in an identifier's value breaking its query
* Fixed the editor reading the schema two or three times when a configuration is set, and showing the tables of the catalog or connection just left until the new ones arrived
* core#1146: Search widgets set up their table view once it is in the document instead of 200 ms later; the Helm renderer no longer waits 200 ms before wiring its clicks
* A connection whose schema cannot be loaded (a database the stand cannot reach) is reported as a warning with its error, no longer also as a console error on the first molecule made current

## 1.7.2 (2026-01-21)

* Fix tooltip texts in join editor
* Add explicit references editor in the configurator, with ability to specify schemas in joins

## 1.7.0 (2025-12-30)

* Support multi schema exploration in a single connection
* Map referenced columns to their original tables to extend linked data
* Fix concurent detection errors
* Resolve all correct db information from columns

## 1.6.0 (2025-12-15)

* Support other providers than default PG provider in connection details

## 1.5.0 (2025-12-15)

* Support generic db exploration without registered schema

## 1.4.0 (2025-12-14)

* Refactoring
* Added visual editor for explorer configs

## 1.3.4 (2025-12-11)

* In place registration of regexp detectors for semantic types

## 1.3.0 (2025-10-16)

* Lazy loading of schemas and references

## 1.1.4 (2025-04-10)

* Ability to specify nqName and provider name in connection details

## 1.1.3 (2024-10-03)

Add detection of variable names in search integrated queries

## 1.1.2 (2024-09-27)

Fix value conversion in search results

## 1.1.1 (2024-09-27)

Add pane replacing options

## 1.1.0 (2024-09-27)

Search related features

## 1.0.0 (2024-09-22)

Initial release with features for exploring the database
