from datetime import datetime
from datagrok_api.models.data_connection import DataConnection
from datagrok_api.models.func import Func, FuncParam

@Func.register_subclass("data-query")
class DataQuery(Func):
    """
    Represents a data query function within the Datagrok platform.

    `DataQuery` is a specialized `Func` that encapsulates the concept 
    of querying a data source, such as an external SQL database, using 
    a provided connection. This class combines the platform's general 
    function capabilities with specific query execution details.

    Attributes
    ----------
    query : str
        The query text to be executed. Typically SQL, but may vary 
        depending on the data source type. Can contains special annotations 
        that will be parsed by platform into the func parameters
    connection : DataConnection
        The connection definition for the data source on which the 
        query will be executed.
    **kwargs
        Additional parameters passed to the base `Func` class, such as:
            - id : Optional[str]
            - name : Optional[str]
            - friendly_name : Optional[str]
            - description : Optional[str]
            - created_on : Optional[datetime]
            - updated_on : Optional[datetime]
            - params : Optional[List[FuncParam]]
            - namespace : Optional[str]
            - tags : Optional[List[str]]
            - options : Optional[Dict[str, Any]]
    query : str
        The raw query text.
    connection : DataConnection
        Connection settings for the target data source.

    Examples
    --------
    Creating and serializing a `DataQuery`:

    >>> conn = DataConnection(name="Postgres", server="db.example.com", port=5432)
    >>> dq = DataQuery(query="SELECT * FROM customers", connection=conn, name="Customer Query")
    >>> dq.to_dict()
    ... {
    ...     "id": None,
    ...     "name": "Customer Query",
    ...     "friendlyName": None,
    ...     "createdOn": None,
    ...     "updatedOn": None,
    ...     "source": "data-query",
    ...     "description": None,
    ...     "params": [],
    ...     "namespace": None,
    ...     "tags": [],
    ...     "options": {},
    ...     "query": "SELECT * FROM customers",
    ...     "connection": {...}
    ... }
    """
    SOURCE = "data-query"

    def __init__(self, query: str, connection: DataConnection, **kwargs):
        kwargs.pop("source", None)
        super().__init__(source=DataQuery.SOURCE, **kwargs)
        self.query = query
        self.connection = connection

    def to_dict(self):
        d = super().to_dict()
        d["query"] = self.query
        d["connection"] = self.connection.to_dict()
        return d
    
    @classmethod
    def _from_dict(cls, data: dict) -> "DataQuery":
        return cls(
            query=data.get("query"),
            connection=DataConnection.from_dict(data["connection"]) if data.get("connection") else None,
            id=data.get("id"),
            name=data.get("name"),
            friendly_name=data.get("friendlyName"),
            description=data.get("description"),
            created_on=datetime.fromisoformat(data["createdOn"]) if data.get("createdOn") else None,
            updated_on=datetime.fromisoformat(data["updatedOn"]) if data.get("updatedOn") else None,
            source=data.get("source"),
            params=[FuncParam.from_dict(d) for d in data.get("params")] if data.get("params") else [],
            namespace=data.get("namespace"),
            tags=data.get('tags', []),
            options=data.get("options")
        )
    