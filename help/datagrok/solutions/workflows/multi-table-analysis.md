---
title: "Multi-table analysis and drill-down"
sidebar_position: 6
description: An end-to-end workflow for analyzing several related tables together in Datagrok, from a quick lookup to a saved master-detail dashboard with drill-down from summaries to rows.
keywords:
  - multiple tables
  - master-detail
  - drill-down
  - link tables
  - join tables
  - aggregate and drill down
  - database explorer
  - multi-table dashboard
mdx:
  format: mdx
---

Real analyses rarely fit in one table. Orders and their line items, compounds
and their measurements, subjects and their visits: the data lives in several
tables that share key columns. This page walks through the ways to work with
such tables together, from a one-off lookup to a saved dashboard where a
click on a summary row reveals the rows behind it.

The examples use the Northwind demo database (tables `orders`,
`order_details`, and `products`), which ships as a demo connection under
**Browse** > **Databases** (for example **PostgresNorthwind** or
**MySQLNorthwind**, depending on your instance), and a results table with one
row per measurement.

## Join or link

Datagrok relates tables in two ways, and the first decision is which one you
need:

* **Join** combines two tables into one wider table. Use it for one-to-one
  relations, for enrichment (adding a few columns from a lookup table), and
  when you want a single flat table to chart. See
  [Join tables](../../../transform/join-tables.md).
* **Link** keeps both tables separate and synchronizes their row state: the
  current row, selection, or filter in one table drives the other. Use it
  for one-to-many relations and for master-detail browsing, where the detail
  side has many rows per master row. See
  [Link tables](../../../transform/link-tables.md).

When both tables come from the same database, prefer joining in the database.
The [Visual Query Editor](../../../access/databases/databases.md#visual-query-editor)
builds the join from foreign keys, and only the result crosses the network.
Join in Datagrok when the sides come from different sources, such as a query
and a spreadsheet.

## Master-detail

The most common multi-table view: pick a row in the master table and see its
related rows in the detail table.

1. Open both tables. For Northwind, expand the connection, right-click
   `orders` and then `order_details`, and select **Get All** for each.
1. On the **Top Menu**, select **Data** > **Link Tables...**
1. Choose `orders` as the first (source) table and `order_details` as the
   second (target) table, pick `orderid` as the key column on both sides,
   and set **Link Type** to `row to filter`.
1. Click **LINK**. On the new link's tab, select **Filter All On No Rows
   Selected**, so that the detail table shows nothing until an order is
   chosen, and click **CLOSE**.

Clicking an order in the `orders` grid now filters `order_details` to its
line items. To see both in one view, go to the `orders` view, add a
[grid](../../../visualize/viewers/grid.md), click its **Gear** icon, and under
**Data** set **Table** to `order_details` and **Row Source** to **Filtered**.
Dock it below the orders grid.

To show the line items inside the orders grid instead, right-click a cell in
`orders` and select **Add** > **Linked Tables** > `order_details`. A column
with the matching rows appears. See
[Data from linked tables](../../../visualize/viewers/grid.md#data-from-linked-tables).

Other link types cover the variations: `selection to filter` shows the
details of several selected masters at once, `filter to filter` propagates
the master's filter panel to the details, and `row to selection` highlights
the matches without hiding anything. For the full list, see
[Link types](../../../transform/link-tables.md#link-types).

## Cascading

Links chain. Add a second link from `order_details` to `products` on
`productid` with the type `filter to filter`, and clicking an order filters
its line items, which in turn filter the products table to the products in
that order. Add a [bar chart](../../../visualize/viewers/bar-chart.md) on
`products` to the orders view (**Table** set to `products`, **Row Source**
set to **Filtered**, split by `categoryid`), and the chart shows the
categories of the products in the selected order.

## Drilling down

Drill-down is the same mechanism applied from a summary to its rows:

* **From an aggregate.** Summarize the measurements with
  [Aggregate rows](../../../transform/aggregate-rows.md), for example count
  and average per compound. Then link the aggregated table to the original
  table on the grouping columns with `row to filter`. Clicking a summary row
  shows the measurements it was computed from. A pivoted table works the
  same way.
* **From a chart.** Any bar chart or pie chart can be the drill-down
  control. Click its **Hamburger** icon and set **On click** to **Filter**,
  and clicking a bar filters the table to that category. See
  [Viewers as filters](../../../visualize/table-view-1.md#viewers-as-filters).
* **From an identifier into the database.** When the
  [Database Explorer](../../../access/databases/databases.md#database-explorer)
  is configured for your database, clicking an identifier such as a compound
  or batch ID anywhere in Datagrok shows the record and everything related
  to it through foreign keys in the **Context Panel**. No detail table needs
  to be loaded and no query written. The related
  [data enrichment](../../../access/databases/databases.md#data-enrichment)
  feature joins columns from the database to the table you are looking at.
* **From a row into a parameterized query.** A query with an input such as
  `compoundId`, run from the **Context Panel** for the current row, returns
  that row's details. See [Retrieve and filter data](retrieve-and-filter.md).

## Multiple views

A table can have several views, each with its own layout, and all views of a
table share the filter and selection. Right-click the table name and select
**Table** > **Add View** to add one. An overview view and a detail view of
the same table stay in step this way. See
[Multiple views](../../../visualize/table-view-1.md#multiple-views).

## Saving as a dashboard

Click **SAVE**. The **Save project** dialog lists all open tables and views.
Links between tables and viewers that point at other tables are saved with
the project. Before you click **OK**, check:

* **Data sync per table.** Each table has its own switch. For the Northwind
  example, keep it on for all three tables, so that the dashboard re-runs the
  three queries on every open. A live query table next to a static reference
  table is a common combination too.
* **Dependencies.** If one table is derived from another (an aggregate, a
  join of two query results), Datagrok records the order and replays it.
  Static tables load first, and then the generation scripts run in the order
  they were recorded.
* **Tables from other projects.** A table opened from another dashboard is
  saved as a **Link** (a read-only reference that follows the original) or a
  **Clone** (an independent copy). Recipients need the **View and use**
  privilege on the original project to see a linked table.

Then share the dashboard. Recipients also need access to the database
connection behind the queries, which shared demo and team connections
usually cover. See [What recipients need](../../concepts/project/dashboard.md#what-recipients-need).

For a live example that needs no database, open the **Table Linking** demo
under **Data Access** in the
[demo app](https://public.datagrok.ai/apps/Tutorials/Demo/Data-Access/Table-Linking).

## Recipes

* **Compare the measurements of five compounds.** Link compounds to
  measurements with `selection to filter`, and select the five rows.
* **See which products a filtered set of orders contains.** Link `orders` to
  `order_details` and `order_details` to `products` with `filter to filter`,
  and use the filter panel on `orders`.
* **Click a bar to see its rows.** Set the bar chart's **On click** to
  **Filter**.
* **Turn a pivot back into rows.** Link the pivoted table to the source on
  the row key with `row to filter`.
* **Look up everything about an ID without loading tables.** Click the ID
  and read the **Context Panel** (requires the Database Explorer to be
  configured).

## Troubleshooting

<details>
<summary>Clicking a master row doesn't filter the detail table</summary>

Check the key columns. Links compare values, so both columns must have the
same type and formatting. Open **Data** > **Link Tables...** to inspect the
links and confirm the link is enabled.

</details>

<details>
<summary>A viewer shows no rows after linking</summary>

Its **Row Source** is set to **Filtered** or **Selected** and the link
produced an empty set. Click a master row that has details, or set **Row
Source** to **All**.

</details>

<details>
<summary>The dashboard opens with one table missing</summary>

The table is dynamic and its source failed, or it is a linked table whose
original project you can't view. See
[Dashboards](../../concepts/project/dashboard.md#troubleshooting).

</details>

See also:

* [Link tables](../../../transform/link-tables.md)
* [Join tables](../../../transform/join-tables.md)
* [Dashboards](../../concepts/project/dashboard.md)
* [Table View](../../../visualize/table-view-1.md)
* [Dashboard lifecycle](dashboard-lifecycle.md)
