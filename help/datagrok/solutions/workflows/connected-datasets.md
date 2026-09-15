---
title: "Work with connected datasets"
sidebar_position: 4
description: An end-to-end workflow for using the databases and file shares already connected to Datagrok, keeping the data current, and building dashboards that update themselves when the source changes.
keywords:
  - connected data
  - internal datasets
  - data sync
  - dynamic dashboard
  - file share
  - data refresh
  - query cache
  - generation script
mdx:
  format: mdx
---

In most organizations, the internal databases and file shares are connected
to Datagrok once, by an administrator or a data owner, and everyone else
consumes them. This page walks through what that looks like from the
consumer's side: finding what is connected, opening data from it, keeping the
data current, and building dashboards that pick up changes in the source
without anyone re-uploading anything.

The examples use a results database (any relational database connected to
Datagrok) and a file share where an instrument drops CSV exports. The demo
Postgres connections and the **Demo** file share under **Browse** work the
same way if you want to try the steps.

## Finding what is connected

Everything connected to Datagrok is in the **Browse** tree:

* **Databases** lists every connection you created or that was shared with
  you, grouped by database type. Expand a connection to see its saved
  queries and, for most databases, its schemas, tables, and columns.
* **Files** lists file shares (S3, Azure, SharePoint, network drives) shared
  with you. Spreadsheets, CSV, and other data files open as tables.
* **Spaces** holds team spaces with the dashboards, queries, and files that
  colleagues published. Each space also has its own file storage.
* **Dashboards** lists dashboards shared with you.
* **My stuff** holds your own entities, favorites, recent items, and what was
  shared directly with you.

To find a dataset by name, tag, or author across all of these, use the search
box at the top of **Browse**. To learn more, see
[Finding available data](../../../access/databases/databases.md#finding-available-data).

Connections, like other entities, are private to their author until shared.
If a database you expect is missing, nobody has shared it with you yet. Ask
the connection owner or your administrator.

## Opening data

There are three ways to get data out of a connected source, and they leave
different traces behind:

* **Run a saved query.** Expand the connection, double-click the query, fill
  in the parameters if it asks for any, and the result opens in a Table
  View. Datagrok records how the table was produced, including the parameter
  values, as the table's _generation script_.
* **Browse a table.** Right-click a table under the connection and select
  **Get TOP 100** to look at it, or **Get All** to load it. This produces a
  generation script too.
* **Open a file.** Double-click a file under **Files** or in a space. The
  table's generation script records the file path. A file opened in another
  way (for example, from a script) has no generation script.

The generation script is what makes everything below possible. It gives the
table its lineage, lets you re-run the source on demand, and lets a
dashboard refresh itself. You can see it in the **Source** pane on the
**Toolbox** of the Table View.

## Keeping data current

A table you opened is a snapshot taken at that moment. It doesn't change
while you look at it, even if the database does. To bring it up to date,
open the **Source** pane on the **Toolbox** and click **Refresh**. For a
parameterized query, the same pane lets you change the parameters before
refreshing.

Two caches can make the refreshed data older than the source:

* **Query cache.** The query author can enable result caching on a query or
  a connection to speed up repeated runs. Until the cache expires or is
  invalidated, a refresh returns the cached result. See
  [Freshness and caching](../../../access/databases/databases.md#freshness-and-caching).
* **File share cache.** A file share can cache file content on a cron
  schedule. A changed file shows up after the cache is flushed, either on
  schedule, manually with the **Clear cache** command on the connection, or
  on every read when **Preflight** is enabled. See
  [Caching](../../../access/files/files.md#caching).

If a refresh shows old data, check these two before suspecting the source.

## Building a self-updating dashboard

A dashboard is a saved Table View: the table plus its viewers, filters, and
formatting. Whether the dashboard stores the data or re-runs the source on
every open is decided by the **Data sync** switch in the **Save project**
dialog:

* **Data sync on** stores only the generation script. Every time the
  dashboard opens, the query re-runs or the file is re-read, and the saved
  layout is applied to the fresh result. The dashboard is always current, but
  it depends on the source being reachable.
* **Data sync off** stores the data as of the save. The dashboard opens
  instantly and works without the source, but it is a snapshot.

For a new dashboard, Data sync is on by default for every table that has a
generation script. To learn more about the trade-offs, see
[Data sync](../../concepts/project/dashboard.md#data-sync).

To build a self-updating dashboard:

1. Open the data from the connected source as described above.
1. Add viewers, filters, and formatting in the Table View.
1. Click **SAVE**, enter a name, check that **Data sync** is on for the
   table, and click **OK**.
1. Share the dashboard. See [Sharing](#sharing) below.

### Example with a database query

A results table in the database grows every day. You save a query
`Results for project X` with a parameter for the project, run it, build a
dashboard with a scatter plot and a filter panel, and save it with Data sync
on. Tomorrow the dashboard opens with today's rows, and the same scatter
plot and filters apply to them. If the dashboard is meant for another
project, the user changes the parameter in the **Source** pane and clicks
**Refresh**, without touching the query.

### Example with a file

A plate reader exports `plate-reader/2026-09/results.csv` to a shared S3
bucket every night. You open the file from **Browse** > **Files**, build a
dashboard, and save it with Data sync on. When tonight's export replaces the
file, anyone who opens the dashboard tomorrow sees tonight's rows. The
exporter must replace the file in place, or write to a stable name such as
`results-latest.csv`, because the generation script references the file by
path. To learn more, see
[Files as a source](../../concepts/project/dashboard.md#files-as-a-source).

## When the source changes

A self-updating dashboard follows the data, but it binds to the source by
name, so changes to the shape of the data matter:

* **New rows** appear on the next open. Nothing to do.
* **New columns** appear in the grid. Existing viewers are unaffected.
* **A column renamed or removed** in the query or the file leaves the viewers
  bound to it showing a placeholder with the missing column name. The
  dashboard still opens. Rebind the viewer or restore the column.
* **A column removed that a recorded transformation step uses** can stop the
  dashboard from opening. Fix the step on the query's **Transformations**
  tab, or restore the column.
* **The query deleted** makes the dashboard fail to open. Recreate the query
  under the same name.
* **The file renamed or moved** breaks the dashboard. Rename it back, or
  rebuild the dashboard on the new path.

Test any change to a shared source by opening the dashboards that depend on
it before the audience does. For the full list and the practices that keep
dashboards robust, see [Changing the source](../../concepts/project/dashboard.md#changing-the-source).

## Sharing

Sharing the dashboard shares its tables and views. When you save with Data
sync on, the query, script, or file connection the table depends on is saved
as part of the dashboard as well, so recipients can re-run it. What recipients
need in addition is access to the **database connection** behind a query,
because the connection is not part of the dashboard. Shared demo and team
connections usually cover this. A private connection makes the dashboard hang
on the loading spinner for the recipient, without an error message.

The simplest arrangement is to keep the dashboard, its queries, and its files
in one [space](../../concepts/project/space.md) and share the space with a
group. Space permissions cascade to everything in it. To learn more, see
[What recipients need](../../concepts/project/dashboard.md#what-recipients-need).

## Recipes

* **See where a dashboard's data comes from.** Open the dashboard and expand
  the **Source** pane on the **Toolbox**.
* **Look at a large table without loading it.** Right-click the table under
  the connection and select **Get TOP 100**.
* **Turn a query result into a dashboard that follows the data.** Run the
  query, build the view, click **SAVE**, keep **Data sync** on.
* **Freeze a dashboard for a report.** Click **SAVE** > **Save a copy** with
  **Data sync** off.
* **Change what a dynamic dashboard shows without editing the query.** Change
  the parameters in the **Source** pane and click **Refresh**.

## Troubleshooting

<details>
<summary>I don't see a connection or a file share that a colleague uses</summary>

It hasn't been shared with you. Ask the owner to share it, or to share the
space it lives in.

</details>

<details>
<summary>The data is stale after a refresh</summary>

Check the query cache and the file share cache first. See
[Keeping data current](#keeping-data-current).

</details>

<details>
<summary>The dashboard never finishes loading for a colleague</summary>

The colleague can't reach the database connection behind a query. Share the
connection. See [Sharing](#sharing).

</details>

See also:

* [Databases](../../../access/databases/databases.md)
* [File shares](../../../access/files/files.md)
* [Dashboards](../../concepts/project/dashboard.md)
* [Retrieve and filter data](retrieve-and-filter.md)
* [Dashboard lifecycle](dashboard-lifecycle.md)
