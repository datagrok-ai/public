---
title: "Layout"
sidebar_position: 2
format: mdx
---

In Datagrok, the visual representation of data (a _layout_) is separate from the data itself
(a [table](../datagrok/concepts/table.md)). This separation lets you save just
the layout and apply it later to a different dataset.

You can save layouts for:
* A [Table View](../datagrok/navigation/views/table-view.md)
* A specific [viewer](viewers/viewers.md)

## Table View layout

A **Table View** layout defines the arrangement of viewers in a **Table View** and their properties. 

To save a **Table View** layout and apply it to a different dataset:
1. Open a table. Add [viewers](viewers/viewers.md), arrange them, and customize the way you want.
2. Save the layout:
   * **In Datagrok**: On the **Top Menu**, click **View > Layout > Save to
   Gallery**. This action saves the layout on a server. 
   * **Locally**: On the **Top Menu**, click **View > Layout > Download**. This action saves the layout to your local drive.
3. Open another table with similar columns.
4. On the **Top Menu**, click **View > Layout > Open Gallery**. This action
   opens a panel with [suggestions](#layout-suggestions). 
5. To apply the layout, find and click the layout you saved earlier. If you
downloaded the layout, drag and drop the layout file into the
view.

To create a copy of your current view, on the  **Top Menu**, click **View >
Layout > Clone**.

![](view-layout.gif)

:::note developers DATAFRAME-SYNCHRONIZED TAGS

Sometimes, you want presentation-related tags (such as color coding) to be
defined on a table level so that all derived **Table Views** will pick it up. To
do this, you need to save these tags as part of the layout. When such a layout
is applied, these tags are re-applied to table/column tags. Such a tag should
start with "%", e.g., "%myColorCoding".

:::

### Layout suggestions

When a **Table View** is open, Datagrok automatically suggests relevant layouts
from the layout gallery, accessible via **Top Menu > View > Layout > Open
Gallery**.

A layout is deemed relevant if its "columns of interest" match a subset of
columns in the current table, determined based on their names, data types,
semantic types, and other metadata.

When saving a layout, the columns of interest are determined as follows:

* For layouts that include non-grid viewers, all data columns used in these
  viewers become the columns of interest.
* For layouts that include only a [grid](viewers/grid.md), all visible
  (non-hidden) columns become the columns of interest.

## Viewer layout

Similarly to a **Table View** layout, you can save and reuse the settings of
individual [viewers](../visualize/viewers/viewers.md). To do this, click a
**Hamburger icon** on the viewer's header and select the command you want under
**General**:

* **Clone**: Add an exact copy of this viewer to your current **Table View**
* **Save to Gallery**: Save to the layouts gallery

To apply the viewer saved earlier, find it in the gallery (**Top Menu** >
**View > Layout > Open Gallery**) and double click it. This adds the viewer to
the **Table View**.

<!--update: Double click will be replaced with the button-->

