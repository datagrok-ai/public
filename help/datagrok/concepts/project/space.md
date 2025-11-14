---
title: "Spaces"
format: mdx
---

Spaces act like folders that contain various [entities](../objects.md), such as tables,
queries, or scripts. 

Spaces are essential for organizing, managing, and sharing data.

## Spaces hierarchy

In Datagrok, there are two types of projects:

* **Root spaces**: Act as the primary _space_ and can include child
  spaces.
* **Child spaces**: Exist under root spaces and are prefixed with the name
  of the root spaces they belong to. For example, the name `Demo:CoffeeCompany` indicates
  that `CoffeeCompany` is a child space under the root space `Demo`. Child spaces inherit [privileges](../../../govern/access-control/access-control.md#permissions) from the root space.

Datagrok automatically creates implicit space projects for
[plugins](../../../develop/how-to/packages/create-package.md) and users:
* **Plugins**: Each plugin version is a child project under the corresponding
  root space.
* **Users**: Unless you move entity to an existing space, any entity you create is
  saved to your personal root space, accessible under **My
  stuff** in the **Browse** view (e.g., `jdoe:MyNewDashboard` or `jdoe:MyNewQuery`). 

## Creating and managing spaces

[Browse](../../navigation/views/browse.md) organizes spaces in a tree
that governs their hierarchy. You can create your own hierarchy under **Namespaces**: 

* **Root spaces**: Right-click **Spaces**, select **Create
Space...**, and name your space in the dialog that opens.
* **Child spaces**: right-click an existing space, select
  **Create Child Space...**, and name your space in the dialog that opens.

You can create as many root and child spaces as you like.

To perform an action on a space, find it in the **Browse** tree and right-click it. This opens the context menu with commands like share, delete, rename, and so on.

## Moving entities between spaces

When you save an entity, you always save it to a space. All newly created
entities are saved to your personal space, visible under **My stuff** on
the **Browse** tree.

If you have the necessary privileges, you can move entities between spaces by
dragging them to a different location in the **Browse** tree. Valid locations are highlighted with a dotted border. Moving entities
impacts their hierarchy, names, and privileges.

![](../../navigation//views/img/namespaces-drag-and-drop.gif)

When moving entities, you have these options:

1. **Clone**:
   * This action creates a copy of the entity in the new space. The entity in the original space remains unaffected.
1. **Move**:
   * The entity is moved to the new space, is automatically renamed, and adopts the permissions of the new space.
   * A view-only copy of the moved entity is created in the original space. This view-only copy is linked to the entity in the new space.
   * Any changes made to the entity in the new space are automatically reflected in the linked copy.
1. **Link**:
   * This action creates a view-only copy of the entity in the new space. This copy is linked to the entity in the original space.
   * Any changes made to the entity in the original space are automatically reflected in the linked copy.

Linked entities are visually distinguished by a **Link** (<FAIcon icon="fa-solid fa-link" size="1x"/>) **icon**. You cannot edit linked entities directly, but you can clone them.

## Removing entities from spaces

To remove an entity from a space simply move it to any other space.

:::danger

Don't to use the **Delete...** command to remove entities from spaces. If you choose the **Delete...** command, it will permanently delete the entity from the server for all users and projects. This action cannot be undone.

:::

## Searching spaces

To find spaces using [smart search](../../navigation/views/browse.md#entity-search), you can use this metadata:

| Field       | Description                            |
|-------------|----------------------------------------|
| name        |                                        |
| description |                                        |
| ID          |                                        |
| createdOn   |                                        |
| updatedOn   |                                        |
| author      | [User](../../../govern/access-control/users-and-groups#users) object |
| starredBy   | [User](../../../govern/access-control/users-and-groups#users) object |
| commentedBy | [User](../../../govern/access-control/users-and-groups#users) object |
| usedBy      | [User](../../../govern/access-control/users-and-groups#users) object |

## See also

* [Browse](../../navigation/views/browse.md)