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

# List all groups that are not personal
api.list_groups(include_personal=False, include_members=True, include_memberships=True)

# Lookup for group by name
developers_group = api.lookup_groups("Developers")[0]
# Request group membership
request = developers_group.request_membership()
# Get all membership requests
requests = developers_group.get_membership_requests()
# Deny or approve requests. You should be admin of the group
if requests[0].from_group.friendly_name == "John Doe" and not requests[0].is_viewed:
	requests[0].approve()
requests[1].deny()

created_group = api.save_group(Group(friendlyName="R&D", description="Bio research group"))

