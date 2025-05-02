# Datagrok python client library

This library can be used for integration with datagrok. It is a python wrapper for public API, that has OpenAPI specification available [here](http://public.datagrok.ai/api/public/api.yaml).

Refer to [Help](http://datagrok.ai/help) for more information.

## Installation

To install package, use [pip](https://pypi.org/project/pip/).

```shell
pip install datagrok-api
```

## Usage

To use API client, import DatagrokClient to your project:

```python
from datagrok_api import DatagrokClient

api = DatagrokClient('your token', 'datagrok url')
```

## License

See [License.md](https://github.com/datagrok-ai/public/blob/master/LICENSE.md).

## Examples

The package uses Pandas for representation of tables and dataframes.

```python
import pandas as pd
from sklearn.datasets import load_iris

# Calls a Datagrok function
# Dataframes, columns and primitive data types are supported
api.call_function('Abs', {
	'x': -3
})

# Downloads file from a given connection
res = api.download_file('system.demofiles', 'demog.csv')
print(res.head())

# Uploads iris dataset as a table
iris_id = api.upload_table('iris', pd.DataFrame(load_iris()['data']))['ID']

# Iris also can be uploaded as a file
# api.upload_file('system.demofiles', 'iris.csv', 'iris.csv')

# Fetches freshly uploaded Iris table 
res = api.download_table(iris_id)
print(res.head())

# Creates dashboard from iris table with layout uploaded from file "iris.layout"
dashboard_id = api.create_dashboard('python-test', iris_id, layout_filename='iris.layout')['ID']

# Shares dashboard with admin
api.share_dashboard(dashboard_id, 'Test')

# Group Management Examples

# List all groups (excluding personal groups)
groups = api.list_groups(include_personal=False)
print("All groups:", [group.friendly_name for group in groups])

# List groups with members and memberships
groups_with_details = api.list_groups(
    include_personal=False,
    include_members=True,
    include_memberships=True
)
for group in groups_with_details:
    print(f"Group: {group.friendly_name}")
    print(f"Members: {[m.friendly_name for m in group.members]}")
    print(f"Memberships: {[m.friendly_name for m in group.memberships]}")

# Search for specific groups
dev_groups = api.lookup_groups("Dev")
print("Development groups:", [group.friendly_name for group in dev_groups])

# Get current user's group
current_user_group = api.get_current_user_group()
print("Current user's group:", current_user_group.friendly_name)

# Create a new group
new_group = Group(friendly_name="Data Science Team", description="Team working on ML projects")
created_group = api.save_group(new_group)
print("Created group:", created_group.friendly_name)

# Get group members
group_members = api.get_group_members(created_group.id)
print(f"Members of {created_group.friendly_name}:", [m.friendly_name for m in group_members])

# Get group memberships
group_memberships = api.get_group_memberships(created_group.id)
print(f"Memberships of {created_group.friendly_name}:", [m.friendly_name for m in group_memberships])

# Request membership in a group
request = api.request_group_membership(created_group.id)
print("Membership request created:", request.id)

# Get membership requests (as group admin)
requests = api.get_group_membership_requests(created_group.id)
for req in requests:
    print(f"Request from {req.from_group.friendly_name}: {'Approved' if req.is_approved else 'Pending'}")

# Approve or deny requests (as group admin)
if requests:
    if requests[0].from_group.friendly_name == "John Doe" and not requests[0].is_viewed:
        api.approve_membership_request(requests[0].id)
    else:
        api.deny_membership_request(requests[0].id)

# Update group information
created_group.description = "Updated description for Data Science Team"
updated_group = api.save_group(created_group, save_relations=True)
print("Updated group:", updated_group.friendly_name, updated_group.description)

# Add Member Examples

# Find John Doe's personal group
john_doe_groups = api.lookup_groups("John Doe")
john_doe_group = next((g for g in john_doe_groups if g.personal), None)
if john_doe_group:
    print(f"Found John Doe's personal group: {john_doe_group.friendly_name}")
    
    # Add John Doe as regular member
    created_group.add_member(john_doe_group)
    print(f"Added {john_doe_group.friendly_name} as regular member to {created_group.friendly_name}")
    
    # Verify the member was added
    updated_group = api.get_group(created_group.id)
    print(f"Members after adding: {[m.friendly_name for m in updated_group.members]}")
    
    # Add John Doe as admin
    created_group.add_member(john_doe_group, is_admin=True)
    print(f"Added {john_doe_group.friendly_name} as admin to {created_group.friendly_name}")
    
    # Verify the admin member was added
    updated_group = api.get_group(created_group.id)
    print(f"Admin members: {[m.friendly_name for m in updated_group.admin_members]}")
    print(f"All members: {[m.friendly_name for m in updated_group.members]}")
else:
    print("John Doe's personal group not found")

```
