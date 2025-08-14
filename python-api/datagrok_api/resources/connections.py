from typing import List, Optional, Union
from datagrok_api.http_client import HttpClient
from datagrok_api.models.data_connection import DataConnection, DataSourceType
from urllib.parse import quote

class ConnectionsClient:
    """
    Client for managing `DataConnection` objects in the Datagrok platform.

    The `ConnectionsClient` provides methods to create, retrieve, update, 
    delete, list, and test data connections via the Datagrok REST API.  
    A data connection represents configuration details needed to connect 
    to an external data source, including authentication parameters.

    This client wraps calls to `/public/v1/connections` endpoints and 
    handles serialization/deserialization of `DataConnection` objects.

    Parameters
    ----------
    client : HttpClient
        An initialized HTTP client used for making API requests.

    Examples
    --------
    Get connection by Grok Name

    >>> from datagrok_api import DatagrokClient
    >>> from datagrok_api.models import DataConnection, Credentials, DatabaseDataSourceType
    >>> grok = DatagrokClient(base_url="https://public.datagrok.ai/api", api_key="Bearer <your-token>")
    >>> conn = grok.connections.get("JohnDoe:MyConnection")
    >>> print(conn.name)
    MyConnection

    Create, save and test DataConnection

    >>> shop_connection = DataConnection(name="TestConn", data_source=DatabaseDataSourceType.Postgres, description="Test conn", server="localhost", 
                                               port=5432, db="shopdb", credentials=Credentials(login="***", password="***"))
    >>> shop_connection_saved = grok.connections.save(shop_connection, save_credentials=True)
    >>> grok.connections.test(shop_connection_saved)

    Get all CLickHouse connections
    >>> click_house_conns = grok.connections.list(smart_filter="dataSource='ClickHouse'")
    """
    def __init__(self, client: HttpClient):
        self.client = client

    def get(self, conn: str) -> DataConnection:
        """Get DataConnection object by id or Grok Name.
        
        Parameters
        ----------
        conn : str
            ID or Grok Name of the connection to retrieve
            
        Returns
        -------
        User
            DataConnection instance with detailed information
        """
        conn = conn.replace(':', '.')
        endpoint = f"/public/v1/connections/{quote(conn)}"
        response = self.client.get(endpoint)
        return DataConnection.from_dict(response.json())

    def save(self, conn: DataConnection, save_credentials: bool = False) -> DataConnection:
        """
        Save a data connection to the backend.

        Parameters
        ----------
        conn : DataConnection
            The data connection object to save.
        save_credentials : bool, optional
            Whether to save the credentials associated with the connection, by default False.

        Returns
        -------
        DataConnection
            The saved data connection returned from the backend, potentially with updated fields.
        """
        conn.ensure_id()
        endpoint = "/public/v1/connections"
        response = self.client.post(endpoint, json=conn.to_dict(), params={'saveCredentials': str(save_credentials).lower()})
        return DataConnection.from_dict(response.json())
    
    def delete(self, conn: Union[DataConnection, str]):
        """
        Delete a data connection by object or ID or Grok Name.

        Parameters
        ----------
        conn : DataConnection or str
            The data connection object or connection ID or Grok Name to delete.

        Raises
        ------
        Exception
            Raises if the HTTP delete request fails.
        """
        if isinstance(conn, DataConnection):
            conn.ensure_id()
            conn = conn.id
        conn = conn.replace(':', '.')      
        endpoint = f"/public/v1/connections/{conn}"
        self.client.delete(endpoint)

    def list(self, smart_filter: Optional[str] = None, 
             data_source: Optional[DataSourceType] = None, tags: Optional[List[str]] = None) -> List[DataConnection]:
        """
        List data connections with optional filters.

        Parameters
        ----------
        smart_filter : str, optional
            Text filter to search connection names or descriptions, by default None.
        data_source : DataSourceType, optional
            Filter connections by data source type, by default None.
        tags : list of str, optional
            Filter connections by associated tags, by default None.

        Returns
        -------
        list of DataConnection
            A list of data connections matching the provided filters.
        """ 
        params = {}
        if smart_filter:
            params["text"] = smart_filter
        if data_source:
            params["dataSource"] = data_source.value
        if tags:
            params["tags"] = ','.join(tags)    
        endpoint = "/public/v1/connections"
        response = self.client.get(endpoint, params=params)
        return [DataConnection.from_dict(data) for data in response.json()]

    def test(self, conn: DataConnection):
        """
        Test if a data connection is available.

        Parameters
        ----------
        conn : DataConnection
            The data connection object to test.

        Raises
        ------
        Exception
            If the connection can't be tested or test fails or the backend returns an error.
        """
        if not DataSourceType.can_be_tested(conn.data_source):
            raise Exception(f"Connection with datasource {conn.data_source.value} can't be tested")

        endpoint = "/public/v1/connections/test"
        response = self.client.post(endpoint, json=conn.to_dict())
        text = response.json()
        if text != "ok":
            raise Exception(f"Connection is not available:\n{text}")
