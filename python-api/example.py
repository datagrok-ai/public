from datagrok_api import DatagrokClient
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
iris_id = api.upload_table('iris', pd.DataFrame(load_iris()['data']))

# Iris also can be uploaded as a file
# api.upload_file('system.demofiles', 'iris.csv', 'iris.csv')

# Fetches freshly uploaded Iris table 
res = api.download_table(iris_id)
print(res.head())

# Creates dashboard from iris table with layout uploaded from file "iris.layout"
dashboard_id = api.create_dashboard('python-test', iris_id, layout_filename='iris.layout')

# Shares dashboard with admin
api.share_dashboard(dashboard_id, 'Test')
