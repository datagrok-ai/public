# Group

Represents a group in the Datagrok system. A group can be a collection of other groups, with hierarchical relationships defined through parent-child relationships. Groups can have members and memberships, and can be personal.

## Attributes

| Name | Type | Description |
| :--- | :--- | :---------- |
| name | Optional[str] | The unique identifier name of the group |
| friendly_name | Optional[str] | A human-readable name for the group |
| description | Optional[str] | A description of the group's purpose |
| created_on | Optional[datetime] | When the group was created |
| updated_on | Optional[datetime] | When the group was last updated |
| personal | bool | Whether this is a personal group |
| hidden | bool | Whether this group is hidden from regular views |
| parents | List[[GroupRelation](Group.md)] | The parent group relationships |
| children | List[[GroupRelation](Group.md)] | The child group relationships |

## Methods

### `__init__()`

Initialize a new Group instance.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| client | Optional[[DatagrokClient](../../datagrok_client/classes/DatagrokClient.md)] | The Datagrok client instance for making API calls **data : dict Dictionary containing group data including id, name, friendlyName, etc. |

### `members()`

Get all non-admin members of this group.

**Returns**

| Type | Description |
| :--- | :---------- |
| List[[Group](Group.md)] | List of non-admin member groups |

### `admin_members()`

Get all admin members of this group.

**Returns**

| Type | Description |
| :--- | :---------- |
| List[[Group](Group.md)] | List of admin member groups |

### `memberships()`

Get all non-admin memberships of this group.

**Returns**

| Type | Description |
| :--- | :---------- |
| List[[Group](Group.md)] | List of non-admin membership groups |

### `admin_memberships()`

Get all admin memberships of this group.

**Returns**

| Type | Description |
| :--- | :---------- |
| List[[Group](Group.md)] | List of admin membership groups |

### `request_membership()`

Request membership in this group.

**Returns**

| Type | Description |
| :--- | :---------- |
| [GroupMembershipRequest](Group.md) | The created membership request |

### `get_membership_requests()`

Get all pending membership requests for this group.

**Returns**

| Type | Description |
| :--- | :---------- |
| List[[GroupMembershipRequest](Group.md)] | List of pending membership requests |

### `add_member()`

Add a member to this group.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| group | [Group](Group.md) | The group to add as a member |
| is_admin | bool, optional | Whether the new member should have admin privileges, by default False |

### `to_dict()`

Convert the group to a dictionary representation.

**Returns**

| Type | Description |
| :--- | :---------- |
| dict | Dictionary containing all group data |

