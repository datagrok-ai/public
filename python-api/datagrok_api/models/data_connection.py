from enum import Enum
from typing import List, Optional, Union
from datetime import datetime
from datagrok_api.models.model import Model, NamedModel


class DataSourceType(str, Enum):
    @classmethod
    def is_file_data_source(cls, source: Union[str, 'DataSourceType']) -> bool:
        if isinstance(source, Enum):
            source = source.value
        return source in FileDataSourceType._value2member_map_

    @classmethod
    def is_database_data_source(cls, source: Union[str, 'DataSourceType']) -> bool:
        if isinstance(source, Enum):
            source = source.value
        return source in DatabaseDataSourceType._value2member_map_

    @classmethod
    def is_system(cls, source: Union[str, 'DataSourceType']) -> bool:
        if isinstance(source, Enum):
            source = source.value
        return source in {
            DatabaseDataSourceType.AWSSM.value,
            DatabaseDataSourceType.PostgresDart.value,
        }

    @classmethod
    def can_be_tested(cls, source: Union[str, 'DataSourceType']) -> bool:
        if isinstance(source, Enum):
            source = source.value
        return source != DatabaseDataSourceType.AWSSM.value
    
    @classmethod
    def from_str(cls, value: str) -> "DataSourceType":
        if DataSourceType.is_file_data_source(value):
            return FileDataSourceType(value)
        elif DataSourceType.is_database_data_source(value):
            return DatabaseDataSourceType(value)
        else:
            raise ValueError(f"Unknown datasource type {value}")


class FileDataSourceType(DataSourceType):
    AzureBlob = 'Azure Blob'
    Dropbox = 'Dropbox'
    GitHub = 'GitHub'
    GoogleCloud = 'GoogleCloud'
    S3 = 'S3'
    EFS = 'Amazon EFS'
    SharePoint = 'SharePoint'
    Files = 'Files'
    CoreWeave = 'CoreWeave'


class DatabaseDataSourceType(DataSourceType):
    Access = 'Access'
    Athena = 'Athena'
    BigQuery = 'BigQuery'
    Cassandra = 'Cassandra'
    DB2 = 'DB2'
    Excel = 'Excel'
    Firebird = 'Firebird'
    HBase = 'HBase'
    Hive = 'Hive'
    Hive2 = 'Hive2'
    MariaDB = 'MariaDB'
    MlFlow = 'MLFlow'
    MongoDB = 'MongoDB'
    MsSql = 'MS SQL'
    MySql = 'MySQL'
    Neo4j = 'Neo4j'
    Oracle = 'Oracle'
    PostgresDart = 'PostgresDart'
    ClickHouse = 'ClickHouse'
    Postgres = 'Postgres'
    Redshift = 'Redshift'
    Snowflake = 'Snowflake'
    SQLite = 'SQLite'
    Socrata = 'Socrata'
    Sparql = 'Sparql'
    Teradata = 'Teradata'
    Twitter = 'Twitter'
    Vertica = 'Vertica'
    Web = 'Web'
    AWSSM = 'AWS'

class Credentials(Model):
    """
    Holds authentication credentials for a data connection.

    Allows dynamic setting and getting of credential fields like login, password,
    access keys, or tokens depending on the data source.

    Attributes
    ----------
    **kwargs : dict
        Arbitrary credential fields passed as key-value pairs.

    Examples
    --------
    Create Credentials with login and password

    >>> creds = Credentials(login="postgres", password="postgres")
    >>> print(creds.login)          # Access login
    >>> creds.password = "new_pass" # Update password dynamically

    Create Credentials with AWS access keys

    >>> aws_creds = Credentials(access_key="AKIA...", secret_key="...")
    >>> print(aws_creds.access_key)     # Access AWS access key
    >>> aws_creds.secret_key = "NEW_SECRET"  # Update secret key dynamically
    """
    _explicit_attrs = {"id", "parameters", "open_parameters"}
    def __init__(self, **kwargs):
        super().__init__(None)
        self.parameters = dict(kwargs)
        self.open_parameters = dict()

    def get(self, key, default=None):
        return self.parameters.get(key, default)
    
    def to_dict(self):
        return {
            "id": self.id,
            "parameters": self.parameters,
            "openParameters": self.open_parameters
        } 

    @classmethod
    def from_dict(cls, data: dict) -> 'Credentials':
        creds = cls(**data.get('parameters', {}))
        creds.id = data.get('id')
        creds.open_parameters = data.get('openParameters', {})
        return creds   

    def __getattr__(self, name):
        if name in self.parameters:
            return self.parameters[name]
        raise AttributeError(f"'Credentials' object has no attribute '{name}'")

    def __getitem__(self, key):
        return self.parameters[key]
    
    def __setitem__(self, key, value):
        self.parameters[key] = value
    
    def __setattr__(self, name, value):
        if name in Credentials._explicit_attrs:
            object.__setattr__(self, name, value)
        else:    
            self.parameters[name] = value

class DataConnection(NamedModel):
    """
    Represents a data connection in Datagrok.

    Attributes
    ----------
        name : Optional[str]
            Internal name of the model, used as an identifier.
        id : Optional[str]
            Unique identifier for the model instance. Can be None for new instances.
        friendly_name : Optional[str]
            Human-readable name for display purposes. Defaults to None.
        created_on : Optional[datetime]
            Timestamp indicating when the model was created. Defaults to None.
        updated_on : Optional[datetime]
            Timestamp indicating when the model was last updated. Defaults to None.
        namespace : Optional[str]
            Optional namespace prefix to distinguish models with the same name. Defaults to None.
        data_source : DataSourceType
            Type of data source this connection uses.
        description : Optional[str]
            Optional textual description of the data connection.
        credentials : Optional[Credentials]
            Authentication credentials used to access the data source.
        tags : Optional[List[str]]
            List of tags associated with the connection for categorization.
        **kwargs : dict
            Additional arbitrary parameters stored in the `parameters` dictionary.

    Examples
    --------
    Create a connection to a Postgres database with login parameters

    >>> conn_pg = DataConnection(
    ...     name="my_postgres_conn",
    ...     data_source=DataSourceType.Postgres,
    ...     server="localhost",
    ...     port=5432,
    ...     db="mydb",
    ...     new Credentials(login="postgres", password="postgres")
    ... )
    >>> print(conn_pg.host)         # Access arbitrary parameter
    >>> conn_pg.password = "new_password"  # Update arbitrary parameter

    Create a connection to an S3 bucket with different parameters

    >>> conn_s3 = DataConnection(
    ...     name="my_s3_conn",
    ...     data_source=DataSourceType.S3,
    ...     bucket="my-bucket",
    ...     region="us-west-2",
    ...     new Credentials(accessKey="AKIA...", secretKey="...")
    ... )
    >>> print(conn_s3.bucket)       # Access bucket name
    >>> conn_s3.region = "us-east-1"  # Update region dynamically
    """
    def __init__(self, 
                 name: str, 
                 data_source: DataSourceType,
                 id: Optional[str] = None,
                 friendly_name: Optional[str] = None,
                 description: Optional[str] = None,
                 credentials: Optional[Credentials] = None,
                 tags: Optional[List[str]] = None,
                 created_on: Optional[datetime] = None,
                 updated_on: Optional[datetime] = None,
                 namespace: Optional[str] = None,
                 **kwargs):
        super().__init__(id=id, name=name, friendly_name=friendly_name, created_on=created_on, updated_on=updated_on, namespace=namespace)
        self.data_source = data_source
        self.description = description
        self.credentials = credentials or Credentials()
        self.tags = tags or []
        self.parameters = dict(kwargs)

    _explicit_attrs = {"id", "name", "data_source", "friendly_name", "description", "credentials", "parameters", "tags", "created_on", "updated_on", "namespace"}


    def get(self, key, default=None):
        return self.parameters.get(key, default)
    
    def to_dict(self):
        return {
            "id": self.id,
            "parameters": self.parameters,
            "name": self.name,
            "friendlyName": self.friendly_name,
            "description": self.description,
            "dataSource": self.data_source.value,
            "credentials": self.credentials.to_dict(),
            "createdOn": self.created_on.isoformat() if self.created_on else None,
            "updatedOn": self.updated_on.isoformat() if self.updated_on else None,
            "tags": self.tags,
            "namespace": self.namespace
        }

    
    @classmethod
    def from_dict(cls, data: dict) -> 'DataConnection':
        credentials = Credentials.from_dict(data.get('credentials', {}))
        conn = cls(
            name=data.get("name") or "",
            data_source=DataSourceType.from_str(data.get("dataSource") if data.get("dataSource") else FileDataSourceType.Files),
            id=data.get("id"),
            friendly_name=data.get("friendlyName"),
            description=data.get("description"),
            credentials=credentials,
            tags=data.get('tags', []),
            namespace=data.get("namespace"),
            **data.get("parameters", {})
        )
        created_str = data.get("createdOn")
        conn.created_on = datetime.fromisoformat(created_str) if created_str else None

        updated_str = data.get("updatedOn")
        conn.updated_on = datetime.fromisoformat(updated_str) if updated_str else None

        return conn    

    def __getattr__(self, name):
        parameters = object.__getattribute__(self, 'parameters')
        if name in parameters:
            return parameters[name]
        raise AttributeError(f"'DataConnection' parameters has no attribute '{name}'")

    def __getitem__(self, key):
        return self.parameters[key]
    
    def __setitem__(self, key, value):
        self.parameters[key] = value
    
    def __setattr__(self, name, value):
        if name in DataConnection._explicit_attrs:
            object.__setattr__(self, name, value)
        else:    
            self.parameters[name] = value
