from datagrok_api.http_client import HttpClient
from datagrok_api.resources.connections import ConnectionsClient
from datagrok_api.resources.files import FilesClient
from datagrok_api.resources.functions import FunctionsClient
from datagrok_api.resources.groups import GroupsClient
from datagrok_api.resources.shares import SharesClient
from datagrok_api.resources.tables import TablesClient
from datagrok_api.resources.users import UsersClient


class DatagrokClient:
    """
    A client for interacting with the Datagrok API.

    This class provides convenient access to various Datagrok resources such as files, tables,
    users, groups, functions, shares, and data connections. It handles authentication and wraps
    lower-level HTTP calls with resource clients.

    Parameters
    ----------
    base_url : str
        The base URL of the Datagrok server (e.g., "https://public.datagrok.ai/api").
    api_key : str
        API key for authentication. It can be obtained from the Datagrok user profile page.
        Should typically start with "Bearer ".

    Attributes
    ----------
    files : FilesClient
        Access to file operations like upload, download, and sync.
    functions : FunctionsClient
        Allows calling and managing Datagrok functions.
    groups : GroupsClient
        Provides access to user and group management operations.
    tables : TablesClient
        Allows reading and uploading Datagrok tables (datasets).
    shares : SharesClient
        Manage entity sharing and permissions.
    users : UsersClient
        Provides information about platform users.
    connections : ConnectionsClient
        Interact with DataConnections (e.g., Postgres, Snowflake, S3, etc.).

    Examples
    --------
    >>> from datagrok_api import DatagrokClient

    >>> grok = DatagrokClient(base_url="https://public.datagrok.ai/api", api_key="Bearer <your-token>")
    
    Download table from server as pd.DataFrame

    >>> tables = grok.tables.download("JohnDoe:MyTable")

    Run a server-side function

    >>> result = grok.functions.call("JohnDoe:MyFunction", {"a": 1, "b": 2})
    >>> print(result)

    Upload a local file

    >>> uploaded_file = grok.files.upload("data.csv")

    Get current user info

    >>> me = grok.users.current()
    >>> print(me.first_name)

    Find and update a group

    >>> group = grok.groups.find("Chemists")[0]
    >>> group.description = "Updated description"
    >>> grok.groups.save(group)

    Create or modify a connection

    >>> from datagrok_api.models import DataConnection, FileDataSourceType, Credentials
    >>> creds = Credentials(accessKey="aws-key", secretKey="aws-secret")
    >>> conn = DataConnection(name="My S3 Bucket", data_source=FileDataSourceType.S3, bucket="my-bucket", region="eu-west-1", credentials=creds)
    >>> conn = grok.connections.save(conn)

    Use as a context manager (auto closes session)

    >>> with DatagrokClient("https://public.datagrok.ai/api", "Bearer <your-token>") as grok:
    ...     print(grok.users.current().name)
    """
    
    def __init__(self, base_url: str, api_key: str):
        self._http = HttpClient(base_url=base_url, api_key=api_key)
        self.files = FilesClient(client=self._http)
        self.functions = FunctionsClient(client=self._http)
        self.groups = GroupsClient(client=self._http)
        self.tables = TablesClient(client=self._http)
        self.shares = SharesClient(client=self._http)
        self.users = UsersClient(client=self._http)
        self.connections = ConnectionsClient(client=self._http)

        
    def close(self):
        self._http.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
