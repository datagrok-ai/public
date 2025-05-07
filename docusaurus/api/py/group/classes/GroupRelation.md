# GroupRelation

Represents a relationship between two groups. This class models the parent-child relationship between groups, including whether the relationship includes admin privileges.

## Attributes

| Name | Type | Description |
| :--- | :--- | :---------- |
| is_admin | bool | Whether this relationship includes admin privileges |
| parent | Optional[[Group](../../group/classes/Group.md)] | The parent group in the relationship |
| child | Optional[[Group](../../group/classes/Group.md)] | The child group in the relationship |

## Methods

### `__init__()`

Initialize a new GroupRelation instance.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| client | [DatagrokClient](../../datagrok_client/classes/DatagrokClient.md) | The Datagrok client instance for making API calls **data : dict Dictionary containing relation data including id, isAdmin, parent, and child |

### `to_dict()`

Convert the group relation to a dictionary representation.

**Returns**

| Type | Description |
| :--- | :---------- |
| dict | Dictionary containing all relation data |

