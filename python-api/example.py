from datagrok_api import DatagrokClient
from datagrok_api.group import Group

import pandas as pd
from sklearn.datasets import load_iris


api = DatagrokClient('YOUR_TOKEN', 'API_URL')

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

