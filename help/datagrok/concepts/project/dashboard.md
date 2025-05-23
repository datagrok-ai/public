---
title: "Dashboards"
sidebar_position: 2
format: mdx
---

Dashboards are [projects](project.md) that contain data (a [dataframe](../table.md)) and the visualizations applied
to it (a [layout](../../../visualize/view-layout.md)). 

<!---
Use dashboards to visually present data in a pre-specified way. In contrast to
[table views](../../navigation/views/table-view.md) that excel at
[exploratory data analysis](../../solutions/workflows/eda.md),
[data wrangling](../../../transform/transform.md), and other table-specific tasks, dashboards trade the ability to quickly
interrogate data in unpredicted ways for delivering the visuals precisely as designed. In particular, here are some
features that are unique to dashboards:

* Visualize data from more than one table at once
* Use pixel-perfect layouts
* Use gadgets

--->

To create a dashboard:

1. Open a table to access the [Table View](../../navigation/views/table-view.md). 
1. In the **Table View**, you can:
   * Add [viewers](../../../visualize/viewers/viewers.md) to visualize your data
   * [Transform data](../../../transform/transform.md) as needed
   * [Add filters](../../navigation/views/table-view.md#select-and-filter)
   * Customize the grid, such as [color-coding grid columns](../../../visualize/viewers/grid.md#color-code-columns)
   * Optionally, add data from [linked tables](../../../transform/link-tables.md):
      * For viewers other than the grid, click the **Gear (<FAIcon icon="fa-solid fa-gear" size="1x" />) icon** and choose the table in the **Context Panel** under **Data** > **Table**.
      * Inside the grid, you can [show data from linked tables in cells](../../../visualize/viewers/grid.md#data-from-linked-tables).
1. [Save](../../navigation/basic-tasks/basic-tasks.md#save-and-share-a-table) your dashboard.
    >Note: When a table is generated by a [function](../../concepts/functions/functions.md) (such as a [database query](../../../access/access.md#data-query)), you can enable the data sync.  With this setting on, the function re-executes each time the project is opened. [Learn more about dynamic dashboards](../../../access/databases/databases.md#creating-dynamic-dashboards-for-query-results).


<!---

## Custom elements

Expand the 'Elements' pane to add gadgets such as a picture, panel, button, etc.

## Custom code

Certain gadgets let you define code that is executed as a reaction to an event, which is typically triggered by a user.
For instance, if you set the Button's 'OnClick' property to `Info("foo")` script, a "foo" message will be shown when a user
clicks on that button.

## Form designer

* Click on an object to select it; its properties appear in [Context Panel](../../navigation/panels/panels.md#context-panel)
* Click-and-drag to select multiple objects at once.

TODO: This feature is not working yet --->


## Resources

YouTube videos:

<div class="help-video-list" style={{display:"flex","flex-wrap":"wrap",}}>

<div class="card" style={{width:"512px",}}>
<iframe src="https://www.youtube.com/embed/TtVjvxMj9Ds?si=8J08Iqbigx2RtR9T" title="YouTube video player" width="512" height="288" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>
  <div class="card-body">
    <h2 class="card-title">Dynamic Dashboards</h2>
    <p class="card-text">Building dynamic dashboards using database queries</p>
  </div>
</div>
</div>
