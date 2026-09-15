---
title: "Dashboard lifecycle"
sidebar_position: 7
description: An end-to-end workflow for publishing a Datagrok dashboard and keeping it alive, covering snapshot versus dynamic data, saving and re-saving shared dashboards, what sharing propagates, and versioning practices.
keywords:
  - dashboard lifecycle
  - data sync
  - sharing a dashboard
  - save a copy
  - personal view customizations
  - dashboard versioning
  - presentation mode
  - retire a dashboard
mdx:
  format: mdx
---

A dashboard goes through the same stages every time: someone builds it, saves
it, decides whether it should follow the data, shares it, changes it, and
eventually replaces or retires it. This page walks through those stages in
order and, at each one, names the choice to make and the practice that
avoids trouble later. For the reference on each feature, see
[Dashboards](../../concepts/project/dashboard.md).

The example is an assay QC dashboard built on a parameterized query and
shared with a team.

## Building

Open the data, either by running a query or by opening a file from
**Browse**, and shape the view: add viewers, a filter panel, color coding,
calculated columns. Until you click **SAVE**, everything lives in your
browser, and closing the tab loses it.

Two practices at this stage pay off later:

* **Keep derived metrics as scalar formulas.** A calculated column such as
  `${signal} / ${background}` survives a missing input with a warning. A
  transformation step or a vector function that loses its input stops the
  dashboard from opening. See
  [Where to put the logic](../../../transform/query-transformations.md#where-to-put-the-logic).
* **Save the layout to the gallery** (**View** > **Layout** > **Save to
  Gallery**) once it looks right. The layout is then a separate entity you
  can reapply to the next dataset or after a bad edit.
* **Name the query with letters, digits, and underscores only.** A dash in
  the query name breaks the generation script, and the dashboard fails to
  load its data on open. Don't rename the table after opening it either, or
  it loses its generation script and can't use Data sync.

## Saving

Click **SAVE**. For a new dashboard, enter a name and, optionally, a
description. It is saved to your personal space under **My stuff**, and you
can move it to a team space later.

The dialog offers two settings that shape the rest of the lifecycle:

* **Data sync**, per table, decides whether the dashboard stores a snapshot
  or re-runs the source on every open. It is on by default for every table
  that has a generation script. See [Snapshot or dynamic](#snapshot-or-dynamic).
* **Presentation mode** hides the toolbox, menus, and panels when the
  dashboard opens. Turn it on for data consumers. They can switch back with
  the **Design mode** link at the top right.

## Snapshot or dynamic

This is the central decision of the lifecycle:

| | Snapshot (Data sync off) | Dynamic (Data sync on) |
|---|---|---|
| Stores | The data, plus the generation script for lineage | The generation script only |
| On open | Shows the data as of the save | Re-runs the query or re-reads the file, then applies the layout |
| Freshness | As old as the last save | Always current |
| Open time | Fast | Depends on the query and the database |
| Without the source | Works | Stays on the loading spinner or fails to open |
| Parameters | Fixed | Editable in the **Source** pane on the **Toolbox** |
| Recipients need | The dashboard | The dashboard, plus access to the database connection |
| Sensitive to source changes | No | Yes. Renamed or removed columns break viewers or the dashboard |

Choose **dynamic** for operational dashboards that must show current data,
such as the assay QC dashboard the team checks every morning. Choose
**snapshot** for reports, review packages, and anything that must stay
reproducible, and for audiences that can't be given access to the source.
One dashboard can mix both, for example a live query table joined to a
static reference table.

To learn more, see [Data sync](../../concepts/project/dashboard.md#data-sync).

## Sharing

Saving doesn't share. Right-click the dashboard in **Browse** and select
**Share...**, enter users or groups, choose the privilege, and click **OK**.

| Privilege | What the recipient can do |
|---|---|
| **View and use** | Open the dashboard, interact with it, download data, save a copy, keep personal view customizations |
| **Full access** | Everything above, plus save the original project, rename, delete, and share it further |

Share with groups rather than with individual users. When the team changes,
you update the group.

Recipients get an in-app notification (or an email, if you entered an email
address), and the dashboard appears under **Browse** > **Dashboards** for
them. A dynamic dashboard opens with the parameters that were in effect when
it was saved, and recipients can change them in the **Source** pane.

### What sharing propagates

Permissions on a project cascade to everything the project contains. When a
dashboard is saved with Data sync on, the query, script, or file connection
each table depends on is saved as part of the dashboard, so sharing the
dashboard lets recipients re-run the query or script and re-read the file.
No separate permission appears on those entities, and they don't show up in
the recipient's Browse tree, but the dashboard works.

What does not travel with the dashboard is the **database connection** behind
a query. Recipients need access to it separately. Shared demo and team
connections usually cover this. If the connection is private to the author,
the recipient's dashboard never finishes loading: the page stays on the
loading spinner with no message.

| Source of a table | Covered by sharing the dashboard | Share separately |
|---|---|---|
| Stored snapshot | Yes | Nothing |
| Database query | The query, yes. The connection, no | The connection, with **View and use** |
| Script | Yes | Nothing |
| File in a file share | Yes, for reading through the dashboard | The folder, if recipients should browse it too |
| Linked table from another project | No | That project, with **View and use** |

The simplest arrangement: keep the dashboard, its queries, and its files in
one [space](../../concepts/project/space.md), and share the space with the
team's group. Space permissions cascade to everything in it, and new queries
saved into the space are shared automatically.

## Re-saving

Once a dashboard is shared, every further **SAVE** asks how to save. The
answer depends on who you are and what you intend:

| You want to... | Choose | Available to |
|---|---|---|
| Publish a change to everyone | **Save original project** | The author and **Full access** recipients |
| Make a variant or a "version" without touching the original | **Save a copy** | Everyone with access |
| Rearrange viewers for yourself only | **Save personal view customizations** | Everyone with access |

**Save a copy** asks, per table, whether to **Clone** it (an independent copy)
or **Link** it (a read-only reference that keeps following the original).
Link when the copy should show the same data as the original, clone when the
copy should be independent.

**Save personal view customizations** stores your layout changes for you
only. The next time you open the dashboard, you see your arrangement, and
everyone else sees the author's. This is the right choice for a recipient
who wants a different set of viewers without forking the dashboard.

## Changing the source

Editing the query behind a dynamic dashboard changes what the dashboard
shows next time. Adding columns and rows is safe, and renaming the query is
safe, because Datagrok rewrites the generation script. Removing or renaming
columns is where dashboards break: viewers bound to a missing column show a
placeholder, and a transformation step that used the column can stop the
dashboard from opening. Deleting the query makes the dashboard fail to open.

Before changing a shared query:

1. Open every dashboard that uses it and note which columns the viewers and
   calculated columns bind to.
1. Make the change, keeping existing column names where possible.
1. Open the dashboards yourself before the audience does.

For the full list of effects, see
[Changing the source](../../concepts/project/dashboard.md#changing-the-source).

## Versioning

Datagrok keeps no history of dashboard versions. Every **Save original
project** overwrites the previous state. The practices that stand in for
version control:

* **Copy before a risky change.** **SAVE** > **Save a copy**, named for its
  purpose, for example `Assay QC (2026-09 rework)`. Once accepted, share the
  copy and retire the original.
* **Snapshot what must not change.** Save review packages and reports with
  Data sync off. A dynamic dashboard can't be reproduced later.
* **Keep the layout in the gallery** and, for a local backup, **View** >
  **Layout** > **Download**.
* **Write the change into the description.** Recipients see it in Browse and
  in the share notification.
* **Prefer copies over edits when the dashboard isn't yours.**

## Retiring

To take a dashboard out of a space without destroying it, move it to another
space, for example your **My stuff**. Moving leaves a view-only linked copy
in the original space, so move that copy too if the team shouldn't see the
dashboard anymore.

To delete a dashboard, right-click it in Browse and select **Delete...**
Deletion is permanent, affects all users, and removes the tables and views
the dashboard owns. Linked tables and the queries it uses stay, because they
belong to other entities. Before deleting, check whether other dashboards
link to its tables (a linked table shows a **Link** icon in Browse).

## Checklist

Before you share a dashboard with a team:

* Data sync is set the way you intend for every table.
* The database connection behind each query is shared with the audience, or
  the dashboard and its sources live in one shared space.
* Column names in the query are stable, the query name has no dashes, and
  derived metrics are scalar formulas.
* The layout is saved to the gallery.
* The description says what the dashboard is for and which parameters
  matter.
* Presentation mode is on if the audience is data consumers.

## Troubleshooting

<details>
<summary>A shared dashboard never finishes loading for a colleague</summary>

The colleague can't access the database connection behind a query. Share the
connection, or move the dashboard and its sources to a shared space.

</details>

<details>
<summary>My layout changes are not visible to others</summary>

You saved personal view customizations. Save again and choose **Save
original project**. You need **Full access** on the dashboard.

</details>

<details>
<summary>The dashboard fails to open after a query change</summary>

The query was deleted, or a transformation step references a column the
query no longer returns. Recreate the query under the same name, or fix the
step on the query's **Transformations** tab.

</details>

See also:

* [Dashboards](../../concepts/project/dashboard.md)
* [Spaces](../../concepts/project/space.md)
* [Layout](../../../visualize/view-layout.md)
* [Work with connected datasets](connected-datasets.md)
* [Multi-table analysis and drill-down](multi-table-analysis.md)
