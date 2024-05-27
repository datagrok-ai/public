# Datagrok python client library

This library can be used for integration with datagrok. It is a python wrapper for public API, that has OpenAPI specification available [here](http://public.datagrok.ai/api/public/api.yaml).

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

```
