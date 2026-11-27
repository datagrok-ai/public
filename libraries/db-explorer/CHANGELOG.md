# db-explorer changelog

## v.next

* Fixed the first card of a page not opening in the context panel when the server spells the connection's namespace in another case than the configuration (a value made before the connection was found kept the configured name)
* Links: the rows that point at a card are grouped in a box headed by what they share with it (`This drug`, `Same target · GROKTGT-000010`), one box for each key the card holds, each read from the column that holds it
* Links: a table is named from the configuration's `entities`, else from its identifier's semantic type and pattern (`DG_BIOLOGICS_DRUG_ID`, `GROKMOL-`: `Drugs (GROKMOL)`), else from its own name without the words most tables of its schema share (`stg_assay_results_v1`: `Assay results`; `sequence_regions` under a sequence: `Regions`); a table that points at the card twice, or at its own table, says how (`Hierarchy · as parent`, `Child compounds`)
* Links: a table's rows are read and counted when its pane opens, for every size; each row is one line, its key and its name, and opens its card in place; the card's own row is tagged `this`; views a curated configuration says nothing about fold under More; the table and column behind a line, and the table's comment, are in its tooltip
* Added `entities` (label, code, plural, description, more, hidden) and the explicit references' `label` and `reverseLabel` to the configuration, and both to the editor (Entity Names, Labels); the editor keeps them through Import JSON and Save
* Fixed a card losing its own links when a column of its row points at its own table (a parent column): the column took the place of the card's key, so a row without a parent showed none and one with a parent showed the parent's
* Fixed the links of a joined table being read with the card's own key instead of the value of the join's column
* Fixed the pane of a table that points at a card twice being renamed after the pane it is in once it had more than 10 rows
* Fixed a value's card cutting its loading indicator while the value loads, as in PowerPack's search for a ChEMBL id
* Fixed an unhandled rejection when a value's associations cannot be read, e.g. a database out of reach whose queries time out: they are left out with a warning, and the card shows its own failure
* Fixed a caller that asked for the schema while it loaded: it waited for ever when the connection was not found, and after a load that failed it still queried the database (and waited for that query to time out too); it now gets no schema, as a later caller does
* Added catalogs: a configuration can name the catalog of its schema (`catalogName`) for a database that has them (Databricks, Snowflake, SQL Server); its tables, joins and references are then read as `catalog.schema.table`, and the editor offers the catalog before the schema and reads the schemas and tables of that catalog alone. A configuration without a catalog reads the connection's own catalog, as before
* Added binding the editor to one connection (`new DBExplorerEditor(onSave, false, connection)`): it then lists no other connection of the server, says which catalog and schema its tables are read from, and Import JSON brings a configuration made elsewhere, such as for the same tables in another catalog, into the schema it is open on
* Fixed the editor's preview cards failing on SQL Server (`limit 1`), and a quote in an identifier's value breaking its query
* Fixed the editor reading the schema two or three times when a configuration is set, and showing the tables of the catalog or connection just left until the new ones arrived
* core#1146: Search widgets set up their table view once it is in the document instead of 200 ms later; the Helm renderer no longer waits 200 ms before wiring its clicks
* Fixed the Helm renderer's drawing vanishing from a card as soon as it appeared (since core#1146): the platform moves a JS input's element into the input's own root a timer tick after the input is made, and the renderer took the element before that move. It now takes it a tick later
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
