---
title: "Projects"
format: mdx
---

Projects act like folders that contain various [entities](../objects.md), such as dataframes,
queries, or scripts. They are essential for organizing, managing, and sharing data.

## Saving entities to projects

When you work on an
entity (like a table or a query), your changes to it aren't saved automatically.
If you close or refresh the browser, all unsaved data will be lost. To save your work, you must manually [upload the modified entity to the
server](../../navigation/basic-tasks/basic-tasks.md#save). 

All Datagrok entities are stored in projects. Subject to your privileges, you may choose between saving an entity to its original project, or saving it to your [personal project](#project-hierarchy).

## Project hierarchy

In Datagrok, there are two types of projects:

* **Root projects**: Act as the primary _namespace_ and can include child
  projects. 
* **Child projects**: Exist under root projects and are prefixed with the name
  of the root project they belong to. For example, the name `Demo:CoffeeCompany` indicates
  that `CoffeeCompany` is a child project under the root project `Demo`.

Datagrok automatically creates root projects for
[plugins](../../../develop/how-to/create-package.md) and users:
* **Plugins**: Each plugin version is a child project under the corresponding
  root project.../../navigation/basic-tasks/img/project-upload-data-sync.png
* **Users**: Unless you choose an existing project, any entity you create is
  saved to your personal root project, accessible under **My
  stuff** in the **Browse** view (e.g., `jdoe:MyNewDashboard` or `jdoe:MyNewQuery`). 
  
The [Browse](../../navigation/views/browse.md) view organizes projects in a tree
that governs their hierarchy. You can create and manage your own project hierarchy under **Namespaces**: 

* To create a custom root project, right-click **Namespaces** and select **Create
Namespace...** This opens a dialog for naming your project.
* To create a child project under an existing project, right-click it and select
  **Create Child Namespace...**. This opens a dialog for naming your project.

You can create as many root and child projects as you like. You can also
manually change the project's type in the [Context Panel](../../navigation/panels/panels.md#context-panel) under
**Namespaces**. Changing a project from child to root moves it to the top
level of the **Namespaces** directory, and updates the names of all entities in
it. 

Subject to your privileges, you can move any entity to your personal root
project or any project within **Namespaces**. To do this, click an object in the
**Browse** tree and start dragging. As you drag, potential directories are
indicated with the dotted border. 

![](../../navigation/views/img/namespaces-drag-and-drop.gif)

Moving entities within **Browse** affects their
hierarchy, names, privileges, and designations as root or child projects. 

When moving entities, you have two options:

* **Move**: Move an entity to a new location. The entity will be renamed, and
  the original location will link to it. The entity will also adopt the
  permissions of the new directory.
* **Link**: Link an entity to its original location. You can't edit linked enities, but you
  can clone them. Linked entities are marked with a link icon.

## Scratchpad

The **Scratchpad** is as a temporary project for entities that haven't been
[saved](../../navigation/basic-tasks/basic-tasks.md#save) in Datagrok or those
that have been modified. It is located at the top of the **Browse** view, just
below the **Top Menu**.

Any entity you open from **Browse** or any table you generate by opening a 
file or running a function like a [database
query](../../../access/databases/databases.md#running-queries), is
automatically added to the Scratchpad.

To permanently save any newly created or modified entities, click the **SAVE**
button located on the **Scratchpad** or at the top of your current [view](../../navigation/views/views.md). This button remains greyed out until additional changes are made.

Entities must be placed in the **Scratchpad** or in another project. Thus, when
you remove an entity from its original project and do not move it to a new one,
it automatically transfers to the **Scratchpad**. If you then remove that entity
from the **Scratchpad**, it closes. To remove entities from the **Scratchpad**
or any open projects, use the **Context Menu**. You can also use drag and drop
to move entities within the **Scratchpad** and between the projects in the **Browse** tree. 

![](scratchpad.gif)

## Searching projects

To find projects using [smart
search](../../navigation/views/browse.md#entity-search), you can use this metadata:

| Field       | Description                            |
|-------------|----------------------------------------|
| name        |                                        |
| description |                                        |
| createdOn   |                                        |
| updatedOn   |                                        |
| author      | [User](../../../govern/user.md) object |
| starredBy   | [User](../../../govern/user.md) object |
| commentedBy | [User](../../../govern/user.md) object |
| usedBy      | [User](../../../govern/user.md) object |

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
