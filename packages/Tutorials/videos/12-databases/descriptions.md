# Databases and queries

Playlist: **Databases and queries**. Feature name of every clip: `How-to · Databases`. Features in
`public/packages/Tutorials/bdd/features/howto/12-databases/`.

All clips use the Biologics connection (Browse > Databases > PostgresDart > Biologics): the demo database the
Biologics package keeps in the server's own Postgres, so no outside database is needed. The tour of the whole
window (01-getting-started/01) only expands this connection, and the Home search box's query sentences and ID
lookups are in 01-getting-started/13 and 14; these clips don't repeat them.

## Plan

One idea per clip, most useful first.

1. Browse a database: connection, schema, tables, a table's preview and its columns, Get All opens it whole.
2. Run a saved query: Run asks for the parameters (choices that come from the database), the result opens with
   the layout saved with the query, the toolbox's Source pane runs it again with another value.
3. Write a SQL query: New Query... on the connection, type SQL, run it, save it; it is listed under the connection.
4. Build a query without SQL: New Visual Query... on a table: group by, aggregate, the result updating as you go.
5. Add a parameter to a query: an `--input:` line and `@name` in the SQL; Run asks for the value.
6. Save viewers with a query: the query editor's Layout tab, viewers from the toolbox, SAVE; running the saved query
   from Browse opens its result with them.
7. See a whole schema at once: Browse on a schema opens the schema view (a box per table with its columns); a box's
   name and columns show in the context panel, its right-click menu is the table's (Get Top 100).

Not made: sharing a query (the Share dialog is the one 06-dashboards-and-sharing/04 films for a dashboard), and a
dashboard kept current by data sync (06-dashboards-and-sharing/06 films the switch; nothing in a clip could show a
reopened dashboard reading new rows, since the demo database does not change).

## 01 - Browse a database (0:49)

Title: Browse a database in Datagrok

Description: How to look through a database from the Browse panel. Under Databases, connections are grouped by
database type; a connection holds its saved queries and its schemas, and a schema lists its tables. Clicking a
table shows a preview of its first rows and its columns in the context panel, and Get All in the table's
right-click menu opens the whole table as a regular table view. Uses the Biologics demo database.

Feature: `public/packages/Tutorials/bdd/features/howto/12-databases/01-browse-a-database.feature`

## 02 - Run a saved query (0:58)

Title: Run a saved database query with parameters in Datagrok

Description: How to run a query someone has already saved. Saved queries are listed under their connection in
Browse; Run in a query's right-click menu asks for its parameters, here an organism and a target picked from lists
that come from the database. The result opens with the layout saved with the query, and the toolbox's Source pane
keeps the parameters, so you can pick another value and press REFRESH to run it again in the same view. Uses the
Antibody Profile query of the Biologics demo database.

Feature: `public/packages/Tutorials/bdd/features/howto/12-databases/02-run-a-saved-query.feature`

## 03 - Write a SQL query (0:43)

Title: Write and save a SQL query in Datagrok

Description: How to write your own SQL against a database connection. New Query... in the connection's right-click
menu in Browse opens the query editor; give the query a name, type the SQL and press the play icon to see the result
under the editor. SAVE stores the query on the server, and it then shows up under its connection like any other
saved query. Uses the Biologics demo database.

Feature: `public/packages/Tutorials/bdd/features/howto/12-databases/03-write-a-sql-query.feature`

## 04 - Build a query without SQL (0:56)

Title: Build a database query without SQL in Datagrok

Description: How to query a database table by picking columns instead of writing SQL. New Visual Query... in a
table's right-click menu in Browse opens the query builder, with rows for Where, Group by, Aggregate and Order by
(Pivot shows once something is aggregated) and the result underneath, which updates after every change. Here the
expression batches table is grouped by expression system, the yield is averaged, and Pivot adds a column for each
chain (HC and LC). Uses the Biologics demo database.

Feature: `public/packages/Tutorials/bdd/features/howto/12-databases/04-build-a-query-without-sql.feature`

## 05 - Add a parameter to a query (0:52)

Title: Add a parameter to a SQL query in Datagrok

Description: How to turn a fixed value in your SQL into a parameter. In the query editor, a line like
--input: string category = "Developability" declares the parameter with its default value, and the SQL refers to it
as @category. When you run the query it asks for the value, starting with the default; here Legacy is typed in and
only the Legacy assay types come back. Uses the Biologics demo database.

Feature: `public/packages/Tutorials/bdd/features/howto/12-databases/05-add-a-parameter-to-a-query.feature`

## 06 - Save viewers with a query (1:04)

Title: Save charts with a database query in Datagrok

Description: How to make a query open its result with the charts you want. In the query editor, the Layout tab shows
the result once you run the query; add viewers from the toolbox, here a box plot and a histogram, and SAVE stores
them with the query. After that, running the query from Browse opens its result with the same viewers. Uses the
expression batches table of the Biologics demo database.

Feature: `public/packages/Tutorials/bdd/features/howto/12-databases/06-save-viewers-with-a-query.feature`

## 07 - See a whole schema at once (0:46)

Title: See all tables of a database schema at once in Datagrok

Description: How to get an overview of a database schema. Browse in a schema's right-click menu (Browse > Databases >
your connection > Schemas) opens the schema view, with a box for every table listing its columns. Clicking a table's
name or a column shows its details in the context panel, and a box's right-click menu has the table's commands; Get
Top 100 opens the first hundred rows as a table. Uses the Biologics demo database.

Feature: `public/packages/Tutorials/bdd/features/howto/12-databases/07-see-a-whole-schema.feature`
