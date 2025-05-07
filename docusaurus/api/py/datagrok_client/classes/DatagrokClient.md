# DatagrokClient

A client for interacting with the Datagrok API. This class provides methods to interact with various Datagrok services including: - Table operations (download/upload) - File operations (download/upload/sync) - Dashboard management - Group management - Function calls The client handles authentication and provides a convenient interface for making API calls to the Datagrok server.

## Attributes

| Name | Type | Description |
| :--- | :--- | :---------- |
| api_key | str | The API key used for authentication. Should start with "Bearer". |
| base_url | str | The base URL of the Datagrok API (e.g., "https://public.datagrok.ai/api") |
| headers | Dict[str, str] | HTTP headers used for authentication in API requests |

## Methods

### `__init__()`

Initialize a new DatagrokClient instance.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| api_key | str | API key for authentication. Can be obtained from user profile page. Should start with "Bearer". |
| base_url | str | Base URL of the Datagrok API (e.g., "https://public.datagrok.ai/api") |

### `download_table()`

Downloads a table from Datagrok.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| name | str | Identifier of a table. Can be accessed from context panel. Several options supppored: UUID, Grok names (e.g. Namespace.Project.Table) |

**Returns**

| Type | Description |
| :--- | :---------- |
| pandas.DataFrame | Dataframe with table data |

### `upload_table()`

Uploads a table to Datagrok.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| name | str | Name for new table. If in grok format, can also specify project and namespace. |
| dataframe | pandas.Dataframe | table data to save |

**Returns**

| Type | Description |
| :--- | :---------- |
| Dict[str, Any] | set of identifiers for the uploaded table |

### `download_file()`

Download a file from Datagrok.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| connector | str |  |
| Identifier of a connector | ID or Grok name. |  |
| path | str | Path to file in a connector |

**Returns**

| Type | Description |
| :--- | :---------- |
| Union[pd.DataFrame, bytes] | If requested file is csv, returns corresponding pd.DataFrame. Otherwise returns bytes with file content. |

### `upload_file()`

Uploads a file to Datagrok.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| connector | str |  |
| Identifier of a connector | ID or Grok name. |  |
| path | str | Path to file in a connector |
| file_path | str | Path to local file that needs to be uploaded |

### `sync_dir()`

Syncs the content of a local directory with a remote directory on a server, updating local files and metadata based on the response. This method compares the current state of the local directory with the remote directory and updates files based on their modified status or whether they have been deleted. It uses a provided connector to connect to a remote service and syncs the data, storing metadata about file versions and deletions in a `.sync_meta.json` file. If a file has been deleted on the server, it will also be deleted locally.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| connector | str |  |
| Identifier of a connector | ID or Grok name. |  |
| remote_path | str | The path on the server connection to sync with. Should be directory |
| dir | str | The local directory to sync. If the directory doesn't exist, it will be created. |
| ext | str, optional | A file extension filter to limit the types of files to sync. If not specified, all files are synced. recursive: bool, optional, default=True Whether to sync the directory recursively, including subdirectories. If False, only files in the root of the directory are synced. sync_meta_file: str, optional, default=".sync_meta.json" The name of the metadata file that stores the eTags and sync information for the files. This file is used to track the current state of the local directory and should be present for syncing to work correctly. Raises Exception If the server response does not contain a multipart body when expected. Notes - If the server returns a `204 No Content` response, it indicates the folder is empty on the remote, so the local directory will be cleared and the sync metadata will be updated accordingly. - If a file is marked as "deleted" on the remote server, it will be removed from the local directory. - The function uses multipart responses for efficient file transfer and parsing of updates. |

### `share_dashboard()`

Shares a dashboard.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| id | str |  |
| Identifier of a project | ID or Grok name. |  |
| groups | str | Comma-separated list of group/user names |
| access | str | Either 'View' or 'Edit' |

### `create_dashboard()`

Create a new dashboard in Datagrok.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| name | str | Name for the new dashboard |
| table_ids | str | Comma-separated list of table IDs to include in the dashboard |
| layout_filename | Optional[str] | Optional path to a JSON file containing the dashboard layout |

**Returns**

| Type | Description |
| :--- | :---------- |
| Dict[str, Any] | Dictionary containing the dashboard identifiers and metadata |

### `call_function()`

Call a Datagrok function with the specified parameters.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| name | str | Name of the function to call |
| invocation_parameters | Dict[str, Any] | Dictionary of parameters to pass to the function |

**Returns**

| Type | Description |
| :--- | :---------- |
| Any | The result of the function call |

### `list_groups()`

List groups from Datagrok with optional filtering and inclusion of related data.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| smart_filter | Optional[str] | Optional smart search filter to apply include_personal : bool, default=False Whether to include personal groups in the results include_members : bool, default=False Whether to include group members in the results include_memberships : bool, default=False Whether to include group memberships in the results |

**Returns**

| Type | Description |
| :--- | :---------- |
| List[[Group](../../group/classes/Group.md)] | List of Group instances matching the criteria |

### `lookup_groups()`

Search for groups by name.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| query | str | Search query to match against group names |

**Returns**

| Type | Description |
| :--- | :---------- |
| List[[Group](../../group/classes/Group.md)] | List of Group instances matching the search query |

### `get_group()`

Get detailed information about a specific group.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| group_id | str | ID of the group to retrieve |

**Returns**

| Type | Description |
| :--- | :---------- |
| [Group](../../group/classes/Group.md) | Group instance with detailed information |

### `get_group_members()`

Get members of a specific group.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| group_id | str | ID of the group to get members for |
| admin | Optional[bool] | If True, returns only admin members. If False, returns only non-admin members. If None, returns all members. |

**Returns**

| Type | Description |
| :--- | :---------- |
| List[[Group](../../group/classes/Group.md)] | List of Group instances representing the members |

### `get_group_memberships()`

Get memberships of a specific group.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| group_id | str | ID of the group to get memberships for |
| admin | Optional[bool] | If True, returns only admin memberships. If False, returns only non-admin memberships. If None, returns all memberships. |

**Returns**

| Type | Description |
| :--- | :---------- |
| List[[Group](../../group/classes/Group.md)] | List of Group instances representing the memberships |

### `request_group_membership()`

Request membership in a group.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| group_id | str | ID of the group to request membership for |

**Returns**

| Type | Description |
| :--- | :---------- |
| [GroupMembershipRequest](../../group/classes/GroupMembershipRequest.md) | The created membership request |

### `get_group_membership_requests()`

Get pending membership requests for a group.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| group_id | str | ID of the group to get membership requests for |

**Returns**

| Type | Description |
| :--- | :---------- |
| List[[GroupMembershipRequest](../../group/classes/GroupMembershipRequest.md)] | List of pending membership requests |

### `approve_membership_request()`

Approve a membership request.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| request_id | str | ID of the membership request to approve |

### `deny_membership_request()`

Deny a membership request.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| request_id | str | ID of the membership request to deny |

### `get_current_user_group()`

Get the current user's group.

**Returns**

| Type | Description |
| :--- | :---------- |
| [Group](../../group/classes/Group.md) | Group instance representing the current user's group |

### `save_group()`

Save or update a group.

**Parameters**

| Name | Type | Description |
| :--- | :--- | :---------- |
| group | [Group](../../group/classes/Group.md) | The group to save or update save_relations : bool, default=False Whether to save group relations (members and memberships) |

**Returns**

| Type | Description |
| :--- | :---------- |
| [Group](../../group/classes/Group.md) | The saved group instance |

