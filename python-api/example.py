"""
Example usage of the DatagrokClient for interacting with the Datagrok platform.

Make sure to install the `datagrok_api` package and replace <your-token> with a valid API token.

API tokens can be generated from the user profile page on your Datagrok server.

Documentation: https://<your-server>/docs/api (replace with real docs link if available)
"""

from datagrok_api import DatagrokClient
from datagrok_api.models import DataConnection, FileDataSourceType, Credentials
import pandas as pd
from sklearn.datasets import load_iris

# Initialize the client
grok = DatagrokClient(
    base_url="https://public.datagrok.ai/api",
    api_key="Bearer <your-token>"  # Replace with your actual API token
)

# Download a Datagrok table as a Pandas DataFrame
table_df = grok.tables.download("JohnDoe:MyTable")
print("Downloaded table:")
print(table_df.head())

# Uploads iris dataset as a table
iris_id = grok.tables.upload('iris', pd.DataFrame(load_iris()['data']))['ID']

# Fetches freshly uploaded Iris table 
res = grok.tables.download(iris_id)
print(res.head())

# Call a Datagrok function
result = grok.functions.call("JohnDoe:MyFunction", {"a": 1, "b": 2})
print("Function result:")
print(result)

# Upload a local file to the Datagrok platform
uploaded_file = grok.files.upload("data.csv")
print(f"Uploaded file: {uploaded_file.name}")

# Retrieve current user information
me = grok.users.current()
print(f"Logged in as: {me.first_name} {me.last_name} ({me.email})")

# Update a group’s description
chemists = grok.groups.find("Chemists")
if chemists:
    group = chemists[0]
    group.description = "Updated description"
    grok.groups.save(group)
    print(f"Updated group '{group.name}'")

# Create or update a connection to an S3 bucket
creds = Credentials(accessKey="aws-key", secretKey="aws-secret")
conn = DataConnection(
    name="My S3 Bucket",
    data_source=FileDataSourceType.S3,
    bucket="my-bucket",
    region="eu-west-1",
    credentials=creds
)
saved_conn = grok.connections.save(conn)
print(f"Saved connection: {saved_conn.name}")

# Recommended: use DatagrokClient as a context manager
with DatagrokClient("https://public.datagrok.ai/api", "Bearer <your-token>") as grok:
    current_user = grok.users.current()
    print(f"Context-managed user: {current_user.name}")
    grok.shares.share
