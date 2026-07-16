---
title: "Domains"
keywords:
  - domain data
  - business objects
  - registration
  - row-level security
  - data catalog
---

:::note

This feature is in Beta

:::

Domains let you keep structured business data — plates, studies, compounds, inventory,
issues — directly in Datagrok, managed like any other platform object. Each domain is a set
of related tables defined by a [plugin](../../develop/how-to/db/domain-schemas.md). Rows can
be searched, browsed, edited, shared with fine-grained permissions, commented on, watched for
changes, and audited — without leaving the platform or setting up a separate database.

## Browsing

To see the registered domains, go to **Browse** > **Platform** > **Domains**. Expand a domain
to see its tables, or click it to see the tables as cards along with a schema diagram showing
how they relate.

Click a table to open it. The table view works like other Datagrok galleries:

* **Search** at the top filters rows as you type. Use plain text, or conditions like
  `status = "open" and title starts "Crash"`.
* **Filters** on the left offer value facets and quick filters (**All**, **Created by me**,
  **Created recently**). Save the current filter and sort as a named preset to reuse it later
  or share it as a URL.
* Switch between **brief**, **card**, and **grid** modes, and change the sort order.
* Large tables show the first 1,000 rows along with the total count — refine the filter to
  narrow down the result. To analyze the full filtered subset, use **Open in Table View** or
  export it to CSV from the ribbon.

Clicking a row shows its details in the
[Context Panel](../../datagrok/navigation/panels/panels.md#context-panel): properties,
sharing, links to related rows, and history. Double-click a row to open its full page — the
row's details, tabs with its related child rows, and its change history. Every row has a
stable URL like `/domains/plates/plate/P-000123` that you can bookmark and share.

To find a row from anywhere, type its identifier into the global search: the general form is
`<domain>.<table>:<key>` (for example, `plates.plate:P-000123`). Plugins can register
friendlier patterns — for example, `GRIT-123` opens the Grit issue directly.

## Creating and editing

Your permissions determine what you can do: the row actions appear only when you have the
corresponding access.

* **Create**: in the table view, click **New...**. The dialog is generated from the table's
  definition: dropdowns for controlled values, pickers for references to other tables and for
  users and groups, and inline validation with the same rules the server enforces. Columns
  you don't have access to are not shown.
* **Edit**: right-click a row and select **Edit...**, or use **Clone...** to create a copy.
* **Edit in grid**: in grid mode, users with Edit access change cells directly. Edits
  accumulate with unsaved-change markers, and **Save** commits them all in one transaction —
  either everything applies or nothing does.

If someone else changed a row while you were editing it, Datagrok shows a conflict dialog:
reload their version or overwrite it.

To change several rows at once, select them and use the bulk **Edit**, **Delete**, or
**Share** actions — each applies as a single all-or-nothing operation.

### Importing and exporting

To load data in bulk, click **Import...** on the table's ribbon. The wizard accepts CSV and
Parquet files, lets you map file columns to table columns, and offers two modes: insert new
rows, or update-or-insert matching them by the table's natural key (for example, a plate
barcode). Before anything is written, a preview reports problems row by row. Choose whether a
single bad row aborts the whole import or the good rows are applied and the rest reported.

To export, use **Open in Table View** (for analysis and visualization) or the CSV download —
both include the full filtered subset, not just the rows on screen.

## Sharing and access

Access is managed with the standard Datagrok
[permissions](../access-control/users-and-groups.md) (View, Edit, Delete, Share). What can be
shared depends on how the plugin defined each table:

* **Table-level**: access is granted on the whole table — anyone with View sees all its rows.
  Typical for reference and lookup tables.
* **Inherited**: rows of a detail table (for example, issues of a project) follow the access
  of the row they belong to. The row's **Effective access** pane explains where its access
  comes from and links to the place where it is managed.
* **Row-level**: individual rows can be shared with users and groups — right-click a row and
  select **Share...**, exactly like sharing any Datagrok entity.

Independently of rows, columns can be restricted: related columns form named groups (for
example, "chemistry" or "procurement" properties), and each group is granted separately.
Users see only the columns their groups are entitled to — everywhere: in grids, dialogs,
exports, and search.

To share a whole domain, right-click it (in the tree or the gallery) and select **Share...** —
grants apply to the domain and all of its tables and column groups at once. To manage access
to a single table, right-click it and select **Share...**, or use the grants editor in its
**Sharing** pane.

## Collaboration

* **Favorites and comments**: star a row or discuss it in its **Chats** pane, like any
  platform object (available for rows that can be individually shared).
* **Watching**: click the bell on a table's ribbon, or **Watch** on a row, to get notified
  when the data changes. Notifications arrive through your usual channels and respect the
  notification preferences in your user profile. Row-level watching is available on tables
  that keep a change history. Bulk imports send a single summary notification to table
  watchers.

## History

Tables with change tracking enabled record every insert, update, and delete together with the
data itself — who changed what, when, and the before/after values. Open a row's **History**
pane to see its timeline. Changes made together (for example, a bulk edit) are grouped, so
related modifications can be traced as one action.

See also:

* [Domain schemas](../../develop/how-to/db/domain-schemas.md) — defining domains in a plugin
* [Users and groups](../access-control/users-and-groups.md)
* [Audit](../audit/audit.md)
* [Sticky meta](sticky-meta.md) — annotating existing objects with custom metadata
