from typing import Any, Dict, List, Optional, Union
from datagrok_api.http_client import HttpClient
from datagrok_api.models.data_connection import DataConnection
from datagrok_api.models.func import Func
from datagrok_api.models.query import DataQuery
from datagrok_api.models.script import Script


class FunctionsClient:
    """
    Client for interacting with Datagrok's Functions API.

    This client provides methods to call, retrieve, list, and create functions
    (including queries and scripts) on a Datagrok instance. It supports direct 
    execution of functions with parameters.

    Parameters
    ----------
    client : HttpClient
        An initialized HTTP client for communicating with the Datagrok API.

    Examples
    --------
    Create Script

    >>> from datagrok_api import DatagrokClient
    >>> from datagrok_api.models import Func, Script, ScriptLanguage, DataQuery
    >>> grok = DatagrokClient(base_url="https://public.datagrok.ai/api", api_key="Bearer <your-token>")
    >>> script = grok.functions.create_script(script="print('Hello from Datagrok!')", name="My Script")
    >>> print(script.script)
    print('Hello from Datagrok!')

    Create DataQuery

    >>> conn = grok.connections.get("JohnDoe:Chembl")
    >>> query = grok.functions.create_query(conn=conn, query="select id from molecules limit 100", name="Molecules")

    Execute query

    >>> res = grok.functions.call(query)

    Get all scripts

    >>> r_scripts = grok.functions.list_scripts()
    """
    def __init__(self, client: HttpClient):
        self.client = client

    def call(self, func: Union[str, Func], parameters: Optional[Dict[str, Any]] = None) -> Any:
        """Call a Datagrok function with the specified parameters.
        
        Parameters
        ----------
        name : str
            Grok Name of the function to call
        parameters : Dict[str, Any]
            Dictionary of parameters to pass to the function
            
        Returns
        -------
        Any
            The result of the function call
        """

        if isinstance(func, Func):
            func = func.grok_name
        func = func.replace(':', '.')
        endpoint = f"/public/v1/functions/{func}/call"
        response = self.client.post(endpoint, json=parameters if parameters else {})
        return response.json()

    def get(self, id: str) -> Func:
        """Get detailed information about a specific Func.
        
        Parameters
        ----------
        id : str
            ID or Grok Name of the Func to retrieve
            
        Returns
        -------
        Func
            Func instance with detailed information
        """
        id = id.replace(':', '.')
        endpoint = f"/public/v1/functions/{id}"
        response = self.client.get(endpoint)
        return Func.from_dict(response.json())

    def list(self, smart_filter: Optional[str] = None) -> List[Func]:
        """List funcs from Datagrok with optional filtering.
        
        Parameters
        ----------
        smart_filter : Optional[str]
            Optional smart search filter to apply
            
        Returns
        -------
        List[Func]
            List of Func instances matching the criteria
        """
        params = {}
        if smart_filter:
            params["text"] = smart_filter

        endpoint = "/public/v1/functions"
        response = self.client.get(endpoint, params=params)
        return [Func.from_dict(data) for data in response.json()]
    
    def list_queries(self) -> List[DataQuery]:
        return self.list(smart_filter=f"source=\"{DataQuery.SOURCE}\"")
    
    def list_scripts(self) -> List[Script]:
        return self.list(smart_filter=f"source=\"{Script.SOURCE}\"")
    
    def create_query(self, connection: DataConnection, query: str, name: Optional[str]=None) -> DataQuery:
        return self._save(DataQuery(query=query, connection=connection, name=name))
    
    def create_script(self, script: str, name: Optional[str]=None) -> Script:
        return self._save(Script(script=script, name=name))
    
    def _save(self, func: Func) -> Func:
        func.ensure_id()
        endpoint = f"/public/v1/functions"
        response = self.client.post(endpoint, json=func.to_dict())
        return Func.from_dict(response.json())
    