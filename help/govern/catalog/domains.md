---
title: "Domains"
keywords:
  - domain data
  - business objects
  - registration
  - row-level security
  - data catalog
  - filtering
---

:::note

This feature is in Beta

:::

Domains let you keep structured business data — plates, studies, compounds, inventory,
issues — directly in Datagrok, managed like any other platform object. Each domain is a set
of related tables defined by a [plugin](../../develop/how-to/db/domain-schemas.md). Rows can
be searched, browsed, filtered, edited, shared with fine-grained permissions, commented on,
watched for changes, and audited — without leaving the platform or setting up a separate
database.

## Browsing

To see the registered domains, go to **Browse** > **Platform** > **Domains**. Expand a domain
to see its tables, or click it to see the tables as cards along with a schema diagram showing
how they relate.

Click a table to open it. The table view works like other Datagrok galleries:

* **Search** at the top filters rows as you type. Use plain text, or conditions like
  `status = "open" and title starts "Crash"`.
* **Filters** on the left narrow the rows down by column values — see
  [Filtering](#filtering).
* Switch between **brief**, **card**, and **grid** modes, and change the sort order.
* Large tables show the first 1,000 rows along with the total count — refine the filter to
  narrow down the result. To analyze the full filtered subset, use **Open in Table View**
  (see [Queries and dashboards](#queries-and-dashboards)) or export it to CSV from the ribbon.

Clicking a row shows its details in the
[Context Panel](../../datagrok/navigation/panels/panels.md#context-panel): properties,
sharing, links to related rows, and history. Double-click a row to open its full page — the
row's details, tabs with its related child rows, and its change history. Every row has a
stable URL like `/domains/plates/plate/P-000123` that you can bookmark and share.

To find a row from anywhere, type its identifier into the global search: the general form is
`<domain>.<table>:<key>` (for example, `plates.plate:P-000123`). Plugins can register
friendlier patterns — for example, `GRIT-123` opens the Grit issue directly.

## Filtering

To open the filter panel, click the **filter** icon next to the search bar. Datagrok builds
the filters from the table's columns automatically (or shows the set the plugin declared):
checkboxes with live row counts for categorical columns, histograms with editable min and max
inputs for numeric columns, date ranges, and a contains input for free-text columns.

Filters always apply to the whole table, not just the rows on screen. Small tables filter
instantly in the browser. Large tables transparently re-query the server as you adjust the
filters, so the counts and histograms reflect all matching rows, not only the loaded ones.
Counts respect your row-level access: two users filtering the same table can see different
counts when they see different rows. To analyze or chart the filtered subset with the full
power of Datagrok, open it as a regular table — see
[Queries and dashboards](#queries-and-dashboards).

The counter above the table shows how many rows currently match (for example, `124 / 2155`).
When nothing matches, the view says "No rows match the current filter" — click
**Clear filters** to start over.

Working with category filters:

* **Search within a filter**: type in the filter's search box to narrow the displayed list.
  Columns with many distinct values show only the most frequent ones ("top 100 shown — search
  to narrow") — search finds the rest.
* **Only this**: click a category name to keep only its rows.
* **Everything except**: Ctrl-click a checkbox to exclude that value. The exclusion applies
  to the whole table, including values beyond the visible list. Check everything back to
  remove the constraint.
* **Reset**: each filter has its own reset icon (hover over its header) that clears just that
  filter. The panel's reset clears them all.

### Related filters

Click **Add related filter...** to filter by columns of the tables a row refers to — for
example, filter order lines by the category of their product (order line → product →
category), up to three links away. Related filters follow reference columns forward, from a
table to the tables it points at.

### Saved filters

To save the current panel as a named preset, open the **Filter Panel** context menu and
select **Save or Apply > Save...**. Reapply or delete presets from the same menu. Presets are
regular Datagrok entities: share one with users or groups like any other entity, and it
appears in their menu, too.

## Queries and dashboards

To analyze or visualize a filtered subset with grids, viewers, and scripts, open it as a
regular Datagrok table: in the table view, click **Open in Table View** on the ribbon.
Datagrok records the query as a `DomainQuery` function call, so the resulting table knows how
to reproduce itself. Save it in a [project](../../datagrok/concepts/project/dashboard.md) with
**Data sync** enabled, and reopening the project re-runs the query against the current data —
with the access rights of whoever opens it applied at that moment.

When the filter matches more rows than the view has loaded, a dialog offers three choices:

* **CANCEL** — do nothing.
* **FIRST n ONLY** — open exactly the n rows the view has loaded. The recorded query keeps
  `limit = n`, so a dashboard built on it is frozen to the first n rows by design.
* **OPEN ALL M** — open every matching row. The recorded query carries no effective row cap
  (only the 10,000,000-row ceiling all domain queries share), so a dashboard built on it
  stays live: rows added later appear the next time the query runs.

The opened table leads with the table's business columns in their declared order; the
system columns (`id`, `version`, `created_on`, `updated_on`, `author_id`) are hidden
behind a `~` name prefix — show them via the column manager's service-columns toggle.
Note that, being hidden, `id` and `version` are absent from CSV exports of an opened
frame unless you unhide them first.

### The DomainQuery function

`DomainQuery` is a regular [function](../../datagrok/concepts/functions/functions.md): run it
from the console, call it from scripts, or edit it in a table's **Source** pane. It queries
one domain table and returns a dataframe — matching rows (select mode), or grouped totals
when `aggregations` or `groupBy` is set (aggregate mode). Row-level and column-level access
apply on every run.

```js
// Select: filtered, sorted rows ('!' = descending)
DomainQuery("northwind", "orders", filters = ["freight > 30"], orderBy = ["!order_date"])

// Aggregate: group and summarize, following a reference column via a join
DomainQuery("northwind", "products", joins = ["category_id"],
  groupBy = ["category_id.name"], aggregations = ["count", "avg(unit_price) as price"])
```

| Parameter      | Meaning                                                                 |
|----------------|-------------------------------------------------------------------------|
| `schema`       | Domain name                                                             |
| `table`        | Table name                                                              |
| `columns`      | Columns to return (select mode; default — all you can see)              |
| `filters`      | Conditions, combined with AND — see [Filter elements](#filter-elements) |
| `joins`        | Reference columns to follow into the referenced table                   |
| `aggregations` | Measures: `count`, `fn(column)`, or `fn(column) as alias`               |
| `groupBy`      | Grouping columns (aggregate mode)                                       |
| `orderBy`      | Sort columns; prefix with `!` for descending                            |
| `limit`        | Maximum rows to return (default and ceiling — 10,000,000)               |
| `offset`       | Rows to skip (select mode)                                              |

`joins = ["category_id"]` adds the referenced row's columns as `category_id.<name>`, usable
in filters, grouping, and aggregations; in aggregate mode you can also sort by them (a
select-mode `orderBy` accepts own columns only). Measure functions are `count`, `sum`, `avg`,
`min`, and `max`. In aggregate mode `limit` caps the number of groups at 10,000.

### Filter elements

Each `filters` element is one condition, written in any of three forms — the elements are
combined with AND:

* **Filter expression** — the same grammar as the search bar: `freight > 30`,
  `status = "open" and title starts "Crash"`. Always quote text values (`country = "France"`,
  not `country = France`).
* **Condition object** — a JSON condition, like
  `{"property": "region", "operator": "=", "value": "Cote d'Azur"}`. Use this form when the
  value contains an apostrophe or quotation mark (as in `Cote d'Azur`) — such characters
  cannot appear in an expression-grammar value.
* **Condition group** — a JSON list of conditions joined by `"and"`/`"or"`, like
  `[{"property": "freight", "operator": ">=", "value": 20}, "and", {"property": "freight", "operator": "<=", "value": 90}]`.
  This is the form recorded range filters take.

### Aggregation results

Aggregate mode returns one row per group with the requested measures. A few things to know:

* The result is a plain summary table: unlike select mode, its columns carry no domain
  metadata (no choice lists, semantic types, or editors).
* Grouping by a reference column directly (for example, `customer_id`) groups by raw row
  identifiers, so you see UUIDs. To group by something readable, add the reference to
  `joins` and group by a display column: `groupBy = ["customer_id.company_name"]`.
* Boolean grouping columns come back as text values `true` and `false` (empty for unset).

### URL parameters

Dashboards built on `DomainQuery` support
[project URL parameters](../../develop/advanced/url-parameters.md#project-parameters): when
uploading the project, map URL names to the query's parameters — including individual filter
elements, offered as `filters[0]`, `filters[1]`, and so on (aliased as `filters_0` by
default — rename them to friendly names like `region`). Opening
`/p/sales.dashboard?region=France` then re-runs the query with that element replaced, and
sharing links from the open dashboard carries the current values along.

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

To export, use **Open in Table View** (for analysis and visualization — see
[Queries and dashboards](#queries-and-dashboards)) or the CSV download — both include the
full filtered subset, not just the rows on screen.

## Sharing and access

Access is managed with the standard Datagrok
[permissions](../access-control/users-and-groups.md) (View, Edit, Delete, Share). What can be
shared depends on how the plugin defined each table:

* **Table-level**: access is granted on the whole table — anyone with View sees all its rows.
  Typical for reference and lookup tables. A row's read-only **Sharing** pane names the table
  that governs its access.
* **Inherited**: rows of a detail table (for example, issues of a project) follow the access of
  the row they belong to. The row's **Sharing** pane explains where its access comes from and
  links to the place where it is managed.
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
