---
title: "Browse"
keywords:
 - navigation
 - file explorer
 - hierarchical tree
 - data catalog
format: mdx
sidebar_position: 2
---

```mdx-code-block
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';
```

**Browse** is a [view](views.md) that organizes all Datagrok objects (like
files, queries, or projects) in a tree. Similar to Windows File Explorer, it lets you search,
preview, manage, and access anything in Datagrok. 

To open **Browse**, click the **Browse** icon on the **Sidebar**.

![](img/browse-main.gif)

## Controls

Within the **Browse** tree:
* Click an object to preview it.
* Double click to open.
* Right click to see context commands.

Use the up/down (↑↓) keys to navigate and left/right (←→) keys to expand a tree node.

The following commands are available from the **Top Menu**:

* **Home**: Opens the **Home Page**
* **Open local file**: Opens a dialog to upload a file from your local drive
* **Open text**: Opens a [window to import text](../../../access/files/files.mdx#)
* **Refresh view**: Refreshes the **Browse** view

You may see view-specific options on top such as the **Create new connection** button for [databases](../../../access/databases/databases.mdx).

## Tree

**Browse** organizes objects in several ways:

* **My stuff**: Items that we think are relevant to you, such as favorites and recently opened projects.
* **Dashboards**: Dashboards
* **Apps**: Installed apps
* **Namespaces**: Here, you can create
   folders that align with organizational use, such as by project or department. Each folder is a Datagrok [project](link) that can be shared with teams.  
* **Files**: [File shares](../../../access/files/files.mdx)
* **Databases**: Connections, queries, database schemas, tables, and columns
* **Plugins**: All installed [plugins](link)
* **Functions**: Queries, scripts, [OpenAPI](../../../access/open-api.md)
* **Platform**: Users, groups, [connectors](../../../access/databases/connectors/connectors.md), notebooks, predictive models, repositories, dockers, [layouts](../../../visualize/view-layout.md), [sticky meta](../../../govern/catalog/sticky-meta.md)

<!--
* **Platform-centric** (**Platform**): This hierarchy provides quick access
  to the platform's connections, repositories, users and groups, and so on.-->

Subject to your privileges, you can move objects from other directories to your personal directory or
**Namespaces**. To do this, click an object and start dragging. As you drag,
drop zones within the **Browse** tree and **Preview** are highlighted:

* On the **Browse** tree, potential directories are indicated by a dotted border.
* In the **Preview**, drop zones are marked by a gray area.

![](img/namespaces-drag-and-drop.gif)

When moving objects, you have two options:

* **Move**: Move an object to a new directory. The object will be renamed, and
  the original location will link to it. The object will also adopt the
  permissions of the new directory.
* **Link**: Link an object to its original location. You can't edit linked objects but you
  can clone them. Linked object are marked with a link icon.

## Home page

When you first log in, the **Browse** shows a landing view called **Home Page**. This view contains
[widgets](../../../visualize/widgets.md) like **Recent projects** or **Usage**.
The widgets you see depend on the [plugins](link) installed. For example, the
**Recent projects** widget is provided by the [PowerPack
package](https://github.com/datagrok-ai/public/tree/master/packages/PowerPack),
and the **Learn** widget appears with the [Tutorials
package](https://github.com/datagrok-ai/public/tree/master/packages/Tutorials).

Above the main area , there is a search ribbon for [searching anything](link) both within and outside the platform.

<br/>

![](../img/home-page.png)

You can customize the appearance of the **Home Page** as follows:
* Toggle [panels](../panels/panels.md) like **Context Panel** or **Context
  Help**. To do this, use the icons located in the right corner of the **Status
  Bar**.
* Choose which widgets to display, arrange their order, and assign custom names.
  To do this, hover over the widget's top to access its controls. You can close
  a widget or adjust its settings using the [Context
  Panel](../panels/panels.md#context-panel). To manage hidden widgets, go to
  **Sidebar** > **Settings** > **Panels**.

:::note developers

You can customize the **Home Page** programmatically.

:::

<!---

## Previews

### Entity search

Use the free-text input that lets define complex queries. Smart search supports
AND and OR operators and parenthesis, so you can combine filters. If you type
single string - search engine will treat it as filter by name. Tags filtering is
supported: #demo will show entities tagged by #demo tag, also you can combine
tags conditions using AND or OR operators. Every entity has properties, that
could be used for filtering. [See more](../concepts/objects.md).

<details>
<summary>Examples</summary>

Unstructured query; looks for 'biologics' in title and description:

```
Biologics
```

Having #demo tag:

```
# Demo
```

Tagged as either either #demo or #chem:

```
# Demo or #chem
```

Created in the last 7 days:

```
createdOn > -1w
```

Complex conditions:

```
(#demo and #chem) or author = "john@google.com"
starredBy = @current or author = @current
```

Created by recently joined users:

```
author.joined > -5d
```

</details>

## Scratchpad

A special section of the **Browse** view right below the **Top Menu** is called
a **Scratchpad**. Here, you can organize and manage your work in progress. For
example, when you open local files in the platform or open tables, they are
added to the **Scratchpad**. In addition, when you open a table or run query, or job,
or script that produces table - it will be added to scratchpad too. You can remove it using "Remove from project"
context menu. Tables and views must be placed in scratchpad or in another project, so, when you remove table or view
from project - it moves to scratchpad, When you remove table from scratchpad - it closes


To save your work, save the objects in the
**Scratchpad** as a project using **SAVE** button. You can then share it with
others. If you don't make any changes, the **SAVE** button is disabled.

For [Table Views](table-view.md), when you save a project, you save both the underlying data (i.e., a [table](../../concepts/table.md)) and its [Layout]. To learn more, see [Saving](link).

**Scratchpad** is a [project](../../collaborate/project.md). This means you can move tables and other objects within, to and from **Scratchpad** using drag-and-drop, like any other **Browse** directory. To exclude objects from the Scratchpad, use its **Context Menu**.

.

[Scratchpad](scratchpad.md) is a special place to start your own project. You can add data by opening data files, using
drag-and-drop or running queries, then save everything and share with other users. You can exclude entities from project
using context menu, or drag them to another project.

To upload project press _upload_ button in toolbar, or use "Upload" context command.

--->