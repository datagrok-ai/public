# Datagrok python client library

This library can be used for integration with datagrok. It is a python wrapper for public API, that has OpenAPI specification available [here](http://public.datagrok.ai/api/public/api.yaml).

Refer to [Help](http://datagrok.ai/help) for more information.

## Installation

To install package, use [pip](https://pypi.org/project/pip/).

```shell
pip install datagrok-api
```

## Usage

### Initialization
To use the API client, import `DatagrokClient`. It is recommended to use it as a context manager to ensure proper connection handling:

```python
from datagrok_api import DatagrokClient

# Initialize with your server URL and API token
with DatagrokClient("https://public.datagrok.ai/api", "Bearer <your-token>") as grok:
    me = grok.users.current()
    print(f"Logged in as: {me.first_name} {me.last_name}")
```

### Working with Data (Tables & Files)
The package uses Pandas for representation of tables and dataframes.

```python
import pandas as pd
from sklearn.datasets import load_iris

# Upload a DataFrame as a table
iris_df = pd.DataFrame(load_iris()['data'])
iris_id = grok.tables.upload('iris', iris_df)['ID']

# Download a table by ID or Name
df = grok.tables.download(iris_id)
print(df.head())

# Upload a local file
uploaded_file = grok.files.upload("data.csv")
```

### Functions & Connections
You can call server-side functions and manage data connections directly through their respective managers.

```python
from datagrok_api.models import DataConnection, FileDataSourceType, Credentials

# Call a Datagrok function
result = grok.functions.call("JohnDoe:MyFunction", {"a": 1, "b": 2})

# Create and save a connection (e.g., S3)
creds = Credentials(accessKey="aws-key", secretKey="aws-secret")
conn = DataConnection(
    name="My S3 Bucket",
    data_source=FileDataSourceType.S3,
    bucket="my-bucket",
    region="us-east-1",
    credentials=creds
)
saved_conn = grok.connections.save(conn)
```

### Group Management
Manage user groups and metadata.

```python
# Search for a group and update its properties
groups = grok.groups.find("Chemists")
if groups:
    group = groups[0]
    group.description = "Updated description for the team"
    grok.groups.save(group)
```

## License

See [License](https://github.com/datagrok-ai/public/blob/master/LICENSE.md).
