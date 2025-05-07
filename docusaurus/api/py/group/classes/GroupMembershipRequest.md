# GroupMembershipRequest

Represents a request for group membership. This class models a request for a group to join another group, including the approval status and relevant timestamps.

## Attributes

| Name | Type | Description |
| :--- | :--- | :---------- |
| approved | Optional[bool] | Whether the request has been approved (None if pending) |
| created_on | Optional[datetime] | When the request was created |
| updated_on | Optional[datetime] | When the request was last updated |
| resolution_date | Optional[datetime] | When the request was resolved (approved/denied) |
| from_group | Optional[[Group](../../group/classes/Group.md)] | The group requesting membership |
| to | Optional[[Group](../../group/classes/Group.md)] | The group being requested to join |

## Methods

### `__init__()`

Initialize a new GroupMembershipRequest instance.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| client | [DatagrokClient](../../datagrok_client/classes/DatagrokClient.md) | The Datagrok client instance for making API calls **data : dict Dictionary containing request data including id, approved, timestamps, etc. |

### `is_viewed()`

Check if the request has been viewed (approved or denied).

**Returns**

| Type | Description |
| :--- | :---------- |
| bool | True if the request has been viewed, False otherwise |

### `approve()`

Approve this membership request.

### `deny()`

Deny this membership request.

### `to_dict()`

Convert the membership request to a dictionary representation.

**Returns**

| Type | Description |
| :--- | :---------- |
| dict | Dictionary containing all request data |

